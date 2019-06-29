#include <vif.hpp>
#include <vif/astro/ds9.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    std::string model_dir;
    std::string psf;
    std::string image;
    std::string mask;
    std::string out;
    double oversample = 1.0;
    double err = 1.0;
    double xshift = 0.0, yshift = 0.0;
    bool reuse_models = true;
    std::string sim_file;
    std::string suffix;
    bool only_allsub = false;
    bool save_models = true;
    vec1s shift_anchors;
    double max_shift = 0.3;
    double shift_dp = 0.1;
    double shift_fit_radius = 1.5;
    bool use_mcmc = false;
    double ra_center = dnan, dec_center = dnan;
    std::string mcmc_dir = "mcmc";

    read_args(argc, argv, arg_list(
        model_dir, psf, image, mask, out, oversample, name(err, "error"), xshift, yshift,
        reuse_models, sim_file, suffix, only_allsub, shift_anchors, save_models, max_shift,
        shift_fit_radius, shift_dp, use_mcmc, mcmc_dir, ra_center, dec_center
    ));

    if (!suffix.empty() && suffix[0] != '-') suffix = '-'+suffix;

    if (sim_file.empty()) {
        sim_file = file::remove_extension(image)+"-sim"+suffix+".reg";
    }

    if (!file::exists(sim_file)) {
        error("could not read simulation regions from '", sim_file, "'");
        return 1;
    }

    if (!shift_anchors.empty()) {
        xshift = 0; yshift = 0;
        reuse_models = false;
        save_models = false;
    }

    fits::header ihdr;
    vec2d p, img;
    fits::read(psf, p);
    fits::read(image, img, ihdr);
    vec2b bad;
    fits::read(mask, bad);

    bad = bad || !is_finite(img);
    img[where(!is_finite(img))] = 0.0;

    astro::wcs iwcs(ihdr);
    double aspix;
    if (!get_pixel_size(iwcs, aspix)) {
        error("could not get pixel size of image");
        return 1;
    }

    print("pixel size of image: ", aspix, " arcsec");

    double y_center = img.dims[0]/2;
    double x_center = img.dims[1]/2;
    if (is_finite(ra_center)) {
        astro::ad2xy(iwcs, ra_center, dec_center, x_center, y_center);
        x_center -= 1.0; y_center -= 1.0;
    }

    vec3d md;
    model_dir = file::directorize(model_dir);
    vec1s models = model_dir+file::list_files(model_dir, "*.fits");
    inplace_sort(models);

    if (models.empty()) {
        error("could not find any model to fit in ", model_dir);
        return 1;
    }

    if (!shift_anchors.empty()) {
        models = models[where(is_any_of(file::get_basename(models), shift_anchors+".fits"))];
        if (models.empty()) {
            error("no model matches '", shift_anchors, "' for the shift anchor");
            return 1;
        }
    }

    vec2d tp;

    std::string mfile = file::remove_extension(out)+"-models"+suffix+".fits";
    if (reuse_models && file::exists(mfile)) {
        fits::input_image omod(mfile);
        omod.read(md);
    } else {
        for (uint_t i : range(models)) {
            print("generating model ", i+1, " (", models[i], ")");

            vec2d tmd;
            fits::header mhdr;
            fits::read(models[i], tmd, mhdr);
            astro::wcs mwcs(mhdr);

            // Normalize to unit flux
            double ftot = total(tmd);
            if (ftot < 0) {
                error("model ", models[i], " is negative");
                return 1;
            } else if (ftot == 0.0) {
                error("model ", models[i], " is zero");
                return 1;
            }

            tmd /= ftot;

            // Make sure it's big enough
            double maspix;
            if (!get_pixel_size(mwcs, maspix)) {
                error("could not get pixel size of model");
                return 1;
            }

            double factor = 1.0;
            if (tmd.dims[0]*maspix < factor*img.dims[0]*aspix ||
                tmd.dims[1]*maspix < factor*img.dims[1]*aspix) {

                uint_t nx = factor*img.dims[1]*aspix/maspix;
                if (nx % 2 == 0) ++nx;
                uint_t ny = factor*img.dims[0]*aspix/maspix;
                if (ny % 2 == 0) ++ny;

                print(" resizing from ", tmd.dims[0], "x", tmd.dims[1], " to ", ny, "x", nx);

                vec2d omd = tmd;
                tmd = replicate(0.0, ny, nx);
                for (uint_t iy : range(omd.dims[0]))
                for (uint_t ix : range(omd.dims[1])) {
                    tmd.safe(iy + ny/2 - omd.dims[0]/2, ix + nx/2 - omd.dims[1]/2) = omd.safe(iy,ix);
                }

                double crpix1, crpix2;
                fits::getkey(mhdr, "CRPIX1", crpix1);
                fits::getkey(mhdr, "CRPIX2", crpix2);
                fits::setkey(mhdr, "CRPIX1", crpix1 + nx/2 - omd.dims[1]/2);
                fits::setkey(mhdr, "CRPIX2", crpix2 + ny/2 - omd.dims[0]/2);
                fits::setkey(mhdr, "NAXIS1", tmd.dims[1]);
                fits::setkey(mhdr, "NAXIS2", tmd.dims[0]);
                mwcs = astro::wcs(mhdr);
            }

            if (tp.empty()) {
                // Regrid PSF to the model's
                fits::header phdr = ihdr;

                double tra, tdec;
                astro::xy2ad(mwcs, mwcs.dims[1]/2+1, mwcs.dims[0]/2+1, tra, tdec);
                fits::setkey(phdr, "NAXIS1", p.dims[1]);
                fits::setkey(phdr, "NAXIS2", p.dims[0]);
                fits::setkey(phdr, "CRPIX1", p.dims[1]/2 + 1);
                fits::setkey(phdr, "CRPIX2", p.dims[0]/2 + 1);
                fits::setkey(phdr, "CRVAL1", tra);
                fits::setkey(phdr, "CRVAL2", tdec);
                if (oversample > 1.0) {
                    double cd1 = 1.0, cd2 = 1.0;
                    if (fits::getkey(phdr, "CDELT1", cd1) && fits::getkey(phdr, "CDELT2", cd2)) {
                        cd1 /= oversample;
                        cd2 /= oversample;
                        fits::setkey(phdr, "CDELT1", cd1);
                        fits::setkey(phdr, "CDELT2", cd2);
                    } else if (fits::getkey(phdr, "CD1_1", cd1) && fits::getkey(phdr, "CD2_2", cd2)) {
                        cd1 /= oversample;
                        cd2 /= oversample;
                        fits::setkey(phdr, "CD1_1", cd1);
                        fits::setkey(phdr, "CD2_2", cd2);
                    }
                }

                astro::wcs pwcs(phdr);

                print(" regrid PSF to model resolution");
                regrid_interpolate_params rp;
                rp.conserve_flux = true;
                // rp.linearize = true; // buggy?
                rp.method = interpolation_method::cubic;
                tp = regrid_interpolate(p, pwcs, mwcs, rp);
                tp[where(!is_finite(tp))] = 0;
            }

            // Convolve model
            print(" convolve model");
            tmd = convolve2d(tmd, tp);
            if (i == 0) {
                fits::write("tmp.fits", tp);
            }

            // Shift it
            if (xshift != 0.0) {
                double crpix = 0;
                fits::getkey(mhdr, "CRPIX1", crpix);
                fits::setkey(mhdr, "CRPIX1", crpix - xshift/maspix);
            }

            if (yshift != 0.0) {
                double crpix = 0;
                fits::getkey(mhdr, "CRPIX2", crpix);
                fits::setkey(mhdr, "CRPIX2", crpix - yshift/maspix);
            }

            mwcs = astro::wcs(mhdr);

            // Regrid it
            print(" regrid model to image resolution");
            if (maspix < aspix) {
                regrid_drizzle_params rpd;
                rpd.linearize = true;
                tmd = regrid_drizzle(tmd, mwcs, iwcs, rpd);
                tmd[where(!is_finite(tmd))] = 0;
            } else {
                regrid_interpolate_params rp;
                rp.conserve_flux = true;
                // rp.linearize = true; // buggy?
                rp.method = interpolation_method::cubic;
                tmd = regrid_interpolate(tmd, mwcs, iwcs, rp);
                tmd[where(!is_finite(tmd))] = 0;
            }

            if (md.empty()) {
                md = replicate(0.0, models.size()+1, tmd.dims);
                md(0,_,_) = 1.0;
            }

            md(i+1,_,_) = tmd;
        }

        if (save_models) {
            fits::output_image omod(mfile);
            omod.write(md);
            omod.write_header(ihdr);
        }
    }

    if (shift_anchors.empty()) {
        print("fitting");
        vec1u idg = where(bad < 0.5);
        auto res = linfit_pack(img[idg], replicate(err, idg.size()),
            reform(md, md.dims[0], md.dims[1]*md.dims[2])(_,idg));

        print("building residual");
        vec2d residual = img;
        for (uint_t j : range(models)) {
            residual -= res.params[j+1]*md(j+1,_,_);
        }

        residual[where(!is_finite(residual))] = 0;

        vec2f flux_sim;
        vec1f flux_sim_chi2;
        vec1f flux_err_sim;
        if (!sim_file.empty()) {
            print("simulations");
            vec<1,ds9::region> sim;
            ds9::read_regions_physical(sim_file, iwcs, sim);

            if (!use_mcmc) {
                auto batch = linfit_pack_batch(replicate(err, idg.size()),
                    reform(md, md.dims[0], md.dims[1]*md.dims[2])(_,idg));

                flux_sim.resize(res.params.size(), sim.size());
                flux_sim_chi2.resize(sim.size());
                for (uint_t i : range(sim.size())) {
                    vec2b smask(residual.dims);
                    ds9::mask_region(sim[i], smask);
                    vec2d noise = (residual - res.params[0])*vec2d(smask);

                    inplace_translate_integer(noise,
                        round(y_center - sim[i].params[1]),
                        round(x_center - sim[i].params[0])
                    );

                    vec2d timg = img + noise;

                    batch.fit(timg[idg]);

                    flux_sim(_,i) = batch.fr.params;
                    flux_sim_chi2[i] = batch.fr.chi2;
                }
            } else {
                // Read MCMC
                print("reading MCMC models");
                vec1s mcmodels = model_dir+mcmc_dir+"/"+file::list_files(model_dir+mcmc_dir, "*.fits");
                if (mcmodels.empty()) {
                    error("no MCMC model found in ", model_dir+mcmc_dir);
                    return 1;
                }

                vec4f mc;
                vec1u mdid;
                vec1s sid = file::get_basename(models);
                vec1s msid = file::get_basename(mcmodels);
                for (uint_t m : range(mcmodels)) {
                    print("generating models from ", mcmodels[m]);
                    uint_t im = where_first(msid[m] == sid);
                    if (im == npos) {
                        error("MC model ", mcmodels[m], " does not correspond to any normal model");
                        return 1;
                    }

                    mdid.push_back(im+1);

                    vec3f data;
                    fits::input_image iimg(mcmodels[m]);
                    iimg.read(data);
                    fits::header mhdr = iimg.read_header();
                    astro::wcs mwcs(mhdr);

                    // Normalize to unit flux
                    for (uint_t k : range(data.dims[0])) {
                        double ftot = total(data.safe(k,_,_));
                        if (ftot < 0) {
                            error("model ", mcmodels[m], "[", k, "] is negative");
                            return 1;
                        } else if (ftot == 0.0) {
                            error("model ", mcmodels[m], "[", k, "] is zero");
                            return 1;
                        }

                        data.safe(k,_,_) /= ftot;
                    }

                    // Make sure it's big enough
                    double maspix;
                    if (!get_pixel_size(mwcs, maspix)) {
                        error("could not get pixel size of model");
                        return 1;
                    }

                    double factor = 1.0;
                    if (data.dims[1]*maspix < factor*img.dims[0]*aspix ||
                        data.dims[2]*maspix < factor*img.dims[1]*aspix) {

                        uint_t nx = factor*img.dims[1]*aspix/maspix;
                        if (nx % 2 == 0) ++nx;
                        uint_t ny = factor*img.dims[0]*aspix/maspix;
                        if (ny % 2 == 0) ++ny;

                        print(" resizing from ", data.dims[1], "x", data.dims[2], " to ", ny, "x", nx);

                        vec3f odat = data;
                        data = replicate(0.0f, data.dims[0], ny, nx);
                        for (uint_t iy : range(odat.dims[1]))
                        for (uint_t ix : range(odat.dims[2])) {
                            data.safe(_, iy + ny/2 - odat.dims[1]/2, ix + nx/2 - odat.dims[2]/2) =
                                odat.safe(_, iy,ix);
                        }

                        double crpix1, crpix2;
                        fits::getkey(mhdr, "CRPIX1", crpix1);
                        fits::getkey(mhdr, "CRPIX2", crpix2);
                        fits::setkey(mhdr, "CRPIX1", crpix1 + nx/2 - odat.dims[2]/2);
                        fits::setkey(mhdr, "CRPIX2", crpix2 + ny/2 - odat.dims[1]/2);
                        fits::setkey(mhdr, "NAXIS1", data.dims[2]);
                        fits::setkey(mhdr, "NAXIS2", data.dims[1]);
                        mwcs = astro::wcs(mhdr);
                    }

                    // Convolve model
                    print(" convolve models");
                    convolver2d kconv(tp);
                    for (uint_t i : range(data.dims[0])) {
                        vec2d tmd = data.safe(i,_,_);
                        kconv.inplace_convolve(tmd);
                        data.safe(i,_,_) = tmd;
                    }

                    // Shift it
                    if (xshift != 0.0) {
                        double crpix = 0;
                        fits::getkey(mhdr, "CRPIX1", crpix);
                        fits::setkey(mhdr, "CRPIX1", crpix - xshift/maspix);
                    }

                    if (yshift != 0.0) {
                        double crpix = 0;
                        fits::getkey(mhdr, "CRPIX2", crpix);
                        fits::setkey(mhdr, "CRPIX2", crpix - yshift/maspix);
                    }

                    mwcs = astro::wcs(mhdr);

                    // Regrid it
                    print(" regrid models to image resolution");
                    if (maspix <= aspix) {
                        regrid_drizzle_params rpd;
                        rpd.linearize = true;
                        data = regrid_drizzle(data, mwcs, iwcs, rpd);
                        data[where(!is_finite(data))] = 0;
                    } else {
                        error("MCMC does not support resampling to higher resolution");
                        return 1;
                    }

                    if (mc.empty()) {
                        mc = replicate(0.0, mcmodels.size(), data.dims);
                    }

                    mc(m,_,_,_) = data;
                }

                uint_t isim = 0;
                vec3d tmd = md;
                flux_sim.resize(res.params.size(), mc.dims[1]);
                flux_sim_chi2.resize(mc.dims[1]);
                auto pg = progress_start(mc.dims[1]);
                for (uint_t i : range(mc.dims[1])) {
                    // Substitute models
                    tmd.safe(mdid,_,_) = mc.safe(_,i,_,_);
                    if (i == 0) {
                        print(mdid);
                        fits::write("test_model.fits", tmd);
                    }

                    // Add noise
                    vec2b smask(residual.dims);
                    ds9::mask_region(sim[isim], smask);
                    vec2d noise = (residual - res.params[0])*vec2d(smask);

                    inplace_translate_integer(noise,
                        round(y_center - sim[isim].params[1]),
                        round(x_center - sim[isim].params[0])
                    );

                    vec2d timg = img + noise;

                    // Refit
                    auto tres = linfit_pack(timg[idg], replicate(err, idg.size()),
                        reform(tmd, tmd.dims[0], tmd.dims[1]*tmd.dims[2])(_,idg));

                    flux_sim(_,i) = tres.params;
                    flux_sim_chi2[i] = tres.chi2;

                    ++isim;
                    if (isim == sim.size()) isim = 0;

                    progress(pg);
                }
            }

            flux_err_sim = partial_stddev(1, flux_sim);
        } else {
            flux_err_sim = replicate(fnan, res.params.size());
        }

        print("saving");
        if (only_allsub) {
            fits::write(out, residual, ihdr);
        } else {
            fits::output_image oimg(file::remove_extension(out)+"-sub.fits");
            oimg.reach_hdu(1);
            oimg.write(residual);
            oimg.write_header(ihdr);
            oimg.write_keyword("EXTNAME", "RES_ALL");

            for (uint_t i : range(models)) {
                residual = img;
                for (uint_t j : range(models)) {
                    if (j == i) continue;
                    residual -= res.params[j+1]*md(j+1,_,_);
                }

                oimg.reach_hdu(i+2);
                oimg.write(residual);
                oimg.write_header(ihdr);
                oimg.write_keyword("EXTNAME", "RES_"+to_string(i+1));
            }

            fits::write_table(out, "chi2", res.chi2, "flux", res.params[1-_], "flux_err", res.errors[1-_], "flux_err_sim",
                flux_err_sim[1-_], "flux_mc", flux_sim(1-_,_), "chi2_mc", flux_sim_chi2, "model", models);
        }
    } else {
        // Mask pixels too far from the anchor source
        vec1d idmax = mult_ids(bad, max_id(md(1,_,_)));

        for (uint_t iy : range(bad.dims[0]))
        for (uint_t ix : range(bad.dims[1])) {
            double d = sqr(iy - idmax[0]) + sqr(ix - idmax[1]);
            if (d > sqr(shift_fit_radius/aspix)) {
                bad.safe(iy,ix) = 1.0;
            }
        }

        fits::write(file::remove_extension(out)+"-shift-mask.fits", bad, ihdr);

        vec1u idg = where(bad < 0.5);
        auto batch = linfit_pack_batch(replicate(err, idg.size()),
            reform(md, md.dims[0], md.dims[1]*md.dims[2])(_,idg));

        linfit_result res;
        res.chi2 = finf;
        int_t best_x = 0, best_y = 0;
        int_t nshift = ceil((max_shift/aspix)/shift_dp);
        auto pg = progress_start(sqr(2*nshift+1));
        for (int_t iy = -nshift; iy <= nshift; ++iy)
        for (int_t ix = -nshift; ix <= nshift; ++ix) {
            vec2d timg = translate_bicubic(img, iy*shift_dp, ix*shift_dp);
            batch.fit(timg[idg]);
            if (batch.fr.chi2 < res.chi2) {
                res = batch.fr;
                best_x = ix;
                best_y = iy;
            }

            progress(pg);
        }

        img = translate_bicubic(img, best_y*shift_dp, best_x*shift_dp);
        vec2d residual = img;
        for (uint_t j : range(models)) {
            residual -= res.params[j+1]*md(j+1,_,_);
        }
        fits::write(file::remove_extension(out)+"-shift-sub.fits", residual, ihdr);

        print("best shift: ", -best_x*shift_dp*aspix, "  ", -best_y*shift_dp*aspix);
    }

    return 0;
}
