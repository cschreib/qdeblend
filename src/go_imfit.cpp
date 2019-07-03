#include <vif.hpp>
#include <vif/astro/ds9.hpp>

using namespace vif;
using namespace vif::astro;

void write_and_spawn(const std::string& filename, const std::string& str) {
    std::ofstream out(filename);
    out << str << "\n";
    out.close();

    spawn("chmod +x "+filename);
    spawn("./"+filename);
}

int vif_main(int argc, char* argv[]) {
    // std::string fit_image = "wfc3-f160w.fits";
    std::string fit_image = "det_stack.fits";
    std::string seg_image = "det_seg.fits";
    std::string base_out_dir = "imfit/";

    file::mkdir(base_out_dir);

    std::string fit_sub_image = base_out_dir+file::remove_extension(fit_image)+"-sub.fits";

    struct {
        vec1u id, x, y, area;
        vec1f dp, size, flux;
    } incat;

    struct {
        vec1u id;
        vec1f x, y;
        vec1f n, re, ell, pa, flux;
    } outcat;

    fits::read_table("det_cat.fits", ftable(incat.id, incat.x, incat.y, incat.area, incat.flux));

    outcat.id = incat.id;
    outcat.x.resize(incat.id.size());
    outcat.y.resize(incat.id.size());
    outcat.n.resize(incat.id.size());
    outcat.re.resize(incat.id.size());
    outcat.ell.resize(incat.id.size());
    outcat.pa.resize(incat.id.size());
    outcat.flux.resize(incat.id.size());

    // Compute guesses and ranges
    incat.flux /= 0.25*incat.area;
    incat.dp = clamp(sqrt(incat.area/dpi), 2, 4);
    incat.size = 0.25*sqrt(incat.area/dpi);

    vec2u gseg;
    fits::read(seg_image, gseg);
    vec2d gimg;
    fits::header ghdr;
    fits::read(fit_image, gimg, ghdr);
    astro::wcs gw(ghdr);

    double aspix;
    if (!astro::get_pixel_size(gw, aspix)) {
        error("could not determine pixel size");
        return 1;
    }

    vec1i hsizes = max(1.5/aspix, 2*sqrt(incat.area/dpi));

    // From DS9 region files
    std::string base = file::remove_extension(fit_image);
    std::string bad_file = base+"-bad.reg";
    if (!file::exists(bad_file)) {
        bad_file = base+"_bad.reg";
    }
    if (!file::exists(bad_file)) {
        bad_file = "";
    }

    if (!bad_file.empty()) {
        vec2b bad(gimg.dims);
        ds9::mask_regions(bad_file, gw, bad);
        gimg[where(bad)] = dnan;
    }

    double rms = 1.48*mad(gimg[where(gseg == 0)]);

    auto fit_done = [&](uint_t i) {
        return file::exists(base_out_dir+"models/"+to_string(incat.id[i])+".fits");
    };

    file::mkdir(base_out_dir+"models");
    file::mkdir(base_out_dir+"models_conv");

    // Subtract already known models
    bool remove_next = false;
    for (uint_t i : range(incat.id)) {
        if (remove_next) {
            file::remove(base_out_dir+"models/"+to_string(incat.id[i])+".fits");
            file::remove(base_out_dir+"models_conv/"+to_string(incat.id[i])+".fits");
        } else if (!fit_done(i)) {
        //     vec2f cmod;
        //     fits::read(base_out_dir+"models_conv/"+to_string(incat.id[i])+".fits", cmod);
        //     gimg -= cmod;
        // } else {
            remove_next = true;
            break;
        }
    }

    vec1b covered(incat.id.dims);

    // Fit the rest
    for (uint_t i : range(incat.id)) {
        if (covered[i]) continue;

        bool fitted = fit_done(i);

        print("### ", incat.id[i], " (already done: ", fitted, ")");

        std::string out_dir = base_out_dir+to_string(incat.id[i])+"/";
        file::mkdir(out_dir);

        // Get cutout
        vec2u seg;
        vec2b mask;
        int_t hsize = hsizes[i];
        int_t x0 = incat.x[i];
        int_t y0 = incat.y[i];

        if (!fitted) {
            vec1u idi, idr;
            subregion(gseg, {{y0-hsize, x0-hsize, y0+hsize, x0+hsize}}, idi, idr);

            seg = replicate(0, 2*hsize+1, 2*hsize+1);
            vec2f img = replicate(dnan, 2*hsize+1, 2*hsize+1);

            seg[idr] = gseg[idi];
            img[idr] = gimg[idi];

            mask = !is_finite(img);

            // Mask sources that are only partly covered or already substracted
            {
                vec1u sb = seg(0,_);
                append(sb, seg(_,0));
                append(sb, seg(seg.dims[0]-1,_));
                append(sb, seg(_,seg.dims[1]-1));
                vec1u ub = unique_values(sb);
                for (uint_t j : ub) {
                    if (j == 0) continue;
                    vec1u ids = where(seg == j);
                    mask[ids] = true;
                    seg[ids] = 0;
                }

                ub = unique_values(seg);
                for (uint_t j : ub) {
                    if (j == 0) continue;
                    uint_t k = where_first(incat.id == j);
                    if (fit_done(k)) {
                        vec1u ids = where(seg == j);
                        mask[ids] = true;
                        seg[ids] = 0;
                    }
                }
            }

            vec1u idb = where(mask);

            img[idb] = 0;
            vec2f err = sqrt(max(img*0.000446339, 0) + sqr(rms));
            err[idb] = 0;

            // Write data for IMFIT
            fits::output_image oimg(out_dir+"img.fits");
            oimg.write(img);
            oimg.write_header(ghdr);
            double crpix1 = 0, crpix2 = 0;
            if (!fits::getkey(ghdr, "CRPIX1", crpix1) ||
                !fits::getkey(ghdr, "CRPIX2", crpix2)) {
                error("could not read astrometry from base image");
                return 1;
            }
            oimg.write_keyword("CRPIX1", crpix1 - x0 + hsize);
            oimg.write_keyword("CRPIX2", crpix2 - y0 + hsize);

            fits::write(out_dir+"seg.fits", seg);
            // fits::write(out_dir+"error.fits", replicate(rms, img.dims));
            fits::write(out_dir+"error.fits", err);
            fits::write(out_dir+"mask.fits", mask);

            vec2d psf;
            fits::read("det_stack_psf.fits", psf);

            int_t phsize = psf.dims[0]/2;
            if (phsize > 2*hsize) {
                psf = psf((phsize-hsize)-_-(phsize+hsize), (phsize-hsize)-_-(phsize+hsize));
            }

            fits::write(out_dir+"psf.fits", psf);
        } else {
            fits::read(out_dir+"seg.fits", seg);
            fits::read(out_dir+"mask.fits", mask);
        }

        // Do individual fits
        auto write_source_init = [&](std::ofstream& out, uint_t j) {
            double dx = double(incat.x[j])-double(incat.x[i]);
            double dy = double(incat.y[j])-double(incat.y[i]);
            out << "X0 " << dx+hsize+1 << " " << dx+hsize+1 - incat.dp[j] << "," << dx+hsize+1 + incat.dp[j] << " # ID=" << incat.id[j] << "\n";
            out << "Y0 " << dy+hsize+1 << " " << dy+hsize+1 - incat.dp[j] << "," << dy+hsize+1 + incat.dp[j] << "\n";
            out << "FUNCTION Sersic\n";
            out << "PA 90.0\n";
            out << "ell 0.2 0,0.95\n";
            out << "n 1.5 0.1,8.0\n";
            out << "I_0 " << incat.flux[j] << "\n";
            out << "r_e " << incat.size[j] << " 0.01," << incat.size[j]*10 << "\n\n";
        };

        file::mkdir(out_dir+"indiv/");
        vec1u uq = unique_values(seg[where(seg != 0)]);

        vec1f c_x0, c_y0;
        vec1f c_pa, c_ell, c_n, c_ie, c_re;

        for (uint_t nid : uq) {
            std::string sid = to_string(nid);
            print("#### ", sid);

            uint_t j = where_first(incat.id == nid);

            // Fit this source
            if (!fitted) {
                vec2b tmask = mask;
                tmask[where(seg != nid && seg != 0)] = 1;
                fits::write(out_dir+"indiv/mask_"+sid+".fits", tmask);

                std::ofstream out(out_dir+"indiv/imfit_config_"+sid+".dat");

                out << "X0 0 fixed\n";
                out << "Y0 0 fixed\n";
                out << "FUNCTION FlatSky\n";
                out << "I_sky 0 " << -3*rms << "," << 3*rms << "\n\n";

                write_source_init(out, j);
                out.close();

                write_and_spawn(out_dir+"indiv/fit_"+sid+".sh",
                    "/home/cschreib/programming/imfit/imfit \\\n"
                    "-c "+out_dir+"indiv/imfit_config_"+sid+".dat "+out_dir+"img.fits --ftol 1e-6 \\\n"
                    "--mask "+out_dir+"indiv/mask_"+sid+".fits \\\n"
                    "--noise "+out_dir+"error.fits --psf "+out_dir+"psf.fits \\\n"
                    "--save-params "+out_dir+"indiv/imfit_bestfit_"+sid+".dat \\\n"
                    "--save-residual "+out_dir+"indiv/imfit_bestfit_res_"+sid+".fits");
            }

            std::ifstream in(out_dir+"indiv/imfit_bestfit_"+sid+".dat");
            std::string line;
            bool skyfound = false;
            while (std::getline(in, line)) {
                line = trim(line);
                if (line.empty() || line[0] == '#') continue;

                if (line == "FUNCTION FlatSky") {
                    skyfound = true;
                    continue;
                }

                if (!skyfound) continue;

                vec1s spl = split_any_of(line, " \t");
                if (spl.size() != 5) continue;

                float val;
                if (!from_string(spl[1], val)) {
                    error("could not read value for ", spl[0], ": ", spl[1]);
                    return 1;
                }

                if (spl[0] == "X0") {
                    c_x0.push_back(val);
                } else if (spl[0] == "Y0") {
                    c_y0.push_back(val);
                } else if (spl[0] == "PA") {
                    c_pa.push_back(val);
                } else if (spl[0] == "ell") {
                    c_ell.push_back(val);
                } else if (spl[0] == "n") {
                    c_n.push_back(val);
                } else if (spl[0] == "I_e") {
                    c_ie.push_back(val);
                } else if (spl[0] == "r_e") {
                    c_re.push_back(val);
                }
            }
        }

        auto write_source_reuse = [&](std::ofstream& out, uint_t k, uint_t j) {
            out << "X0 " << c_x0[k] << " " << c_x0[k] - incat.dp[j] << "," << c_x0[k] + incat.dp[j] << " # ID=" << incat.id[j] << "\n";
            out << "Y0 " << c_y0[k] << " " << c_y0[k] - incat.dp[j] << "," << c_y0[k] + incat.dp[j] << "\n";
            out << "FUNCTION Sersic\n";
            out << "PA " << c_pa[k] << "\n";
            out << "ell " << c_ell[k] << " 0,0.95\n";
            out << "n " << c_n[k] << " 0.1,8.0\n";
            out << "I_0 " << c_ie[k] << "\n";
            out << "r_e " << c_re[k] << " 0.01," << incat.size[j]*10 << "\n\n";
        };

        if (uq.size() > 1) {
            // Do full fit
            if (!fitted) {
                std::ofstream out(out_dir+"imfit_config.dat");

                out << "X0 0 fixed\n";
                out << "Y0 0 fixed\n";
                out << "FUNCTION FlatSky\n";
                out << "I_sky 0 " << -3*rms << "," << 3*rms << "\n\n";

                for (uint_t k : range(uq)) {
                    uint_t nid = uq[k];
                    uint_t j = where_first(incat.id == nid);
                    write_source_reuse(out, k, j);
                }

                out.close();

                write_and_spawn(out_dir+"fit.sh",
                    "/home/cschreib/programming/imfit/imfit \\\n"
                    "-c "+out_dir+"imfit_config.dat "+out_dir+"img.fits --ftol 1e-6 \\\n"
                    "--mask "+out_dir+"mask.fits \\\n"
                    "--noise "+out_dir+"error.fits --psf "+out_dir+"psf.fits \\\n"
                    "--save-params "+out_dir+"imfit_bestfit.dat \\\n"
                    "--save-residual "+out_dir+"imfit_bestfit_res.fits");
            }

            c_x0.clear();
            c_y0.clear();
            c_pa.clear();
            c_ell.clear();
            c_n.clear();
            c_ie.clear();
            c_re.clear();

            std::ifstream in(out_dir+"imfit_bestfit.dat");
            std::string line;
            bool skyfound = false;
            while (std::getline(in, line)) {
                line = trim(line);
                if (line.empty() || line[0] == '#') continue;

                if (line == "FUNCTION FlatSky") {
                    skyfound = true;
                    continue;
                }

                if (!skyfound) continue;

                vec1s spl = split_any_of(line, " \t");
                if (spl.size() != 5) continue;

                float val;
                if (!from_string(spl[1], val)) {
                    error("could not read value for ", spl[0], ": ", spl[1]);
                    return 1;
                }

                if (spl[0] == "X0") {
                    c_x0.push_back(val);
                } else if (spl[0] == "Y0") {
                    c_y0.push_back(val);
                } else if (spl[0] == "PA") {
                    c_pa.push_back(val);
                } else if (spl[0] == "ell") {
                    c_ell.push_back(val);
                } else if (spl[0] == "n") {
                    c_n.push_back(val);
                } else if (spl[0] == "I_e") {
                    c_ie.push_back(val);
                } else if (spl[0] == "r_e") {
                    c_re.push_back(val);
                }
            }
        }

        for (uint_t k : range(uq)) {
            uint_t nid = uq[k];
            std::string sid = to_string(nid);
            uint_t j = where_first(incat.id == nid);

            std::string ofile = base_out_dir+"models/"+sid+".fits";
            std::string ofile_psf = base_out_dir+"models_conv/"+sid+".fits";

            std::ofstream out(out_dir+"indiv/imfit_bestfit_full_"+sid+".dat");

            c_x0[k] += x0-hsize;
            c_y0[k] += y0-hsize;

            write_source_reuse(out, k, j);
            out.close();

            if (!fitted) {
                write_and_spawn(out_dir+"indiv/make_conv_"+sid+".sh",
                    "/home/cschreib/programming/imfit/makeimage "+
                    out_dir+"indiv/imfit_bestfit_full_"+sid+".dat --refimage "+fit_image+" "
                    "--psf "+out_dir+"psf.fits -o "+ofile_psf+" \n"
                    "fitstool "+fit_image+" cpwcs out="+ofile_psf);

                write_and_spawn(out_dir+"indiv/make_"+sid+".sh",
                    "/home/cschreib/programming/imfit/makeimage "+
                    out_dir+"indiv/imfit_bestfit_full_"+sid+".dat --refimage "+fit_image+" "
                    "-o "+ofile+" \n"
                    "fitstool "+fit_image+" cpwcs out="+ofile);
            }

            vec2f cmod;
            fits::read(ofile_psf, cmod);
            gimg -= cmod;

            covered[j] = true;
        }

        fits::write(fit_sub_image, gimg, ghdr);
    }

    // Read all bestfits
    for (uint_t i : range(incat.id)) {
        std::string out_dir = base_out_dir+to_string(incat.id[i])+"/";
        if (!file::exists(out_dir+"imfit_bestfit.dat")) continue;

        vec1u c_id;
        vec1f c_x0, c_y0;
        vec1f c_pa, c_ell, c_n, c_ie, c_re;

        std::ifstream in(out_dir+"imfit_bestfit.dat");
        std::string line;
        bool skyfound = false;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;

            if (line == "FUNCTION FlatSky") {
                skyfound = true;
                continue;
            }

            if (!skyfound) continue;

            vec1s spl = split_any_of(line, " \t");
            if (spl.size() != 5) continue;

            float val;
            if (!from_string(spl[1], val)) {
                error("could not read value for ", spl[0], ": ", spl[1]);
                return 1;
            }

            if (spl[0] == "X0") {
                c_x0.push_back(val);
            } else if (spl[0] == "Y0") {
                c_y0.push_back(val);
            } else if (spl[0] == "PA") {
                c_pa.push_back(val);
            } else if (spl[0] == "ell") {
                c_ell.push_back(val);
            } else if (spl[0] == "n") {
                c_n.push_back(val);
            } else if (spl[0] == "I_e") {
                c_ie.push_back(val);
            } else if (spl[0] == "r_e") {
                c_re.push_back(val);
            }
        }

        // Read IDs...
        in.close();
        in.open(out_dir+"imfit_config.dat");
        skyfound = false;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty() || line[0] == '#') continue;

            if (line == "FUNCTION FlatSky") {
                skyfound = true;
                continue;
            }

            if (!skyfound) continue;

            vec1s spl = split_any_of(line, " \t");
            // print(spl);
            if (spl.size() != 5) continue;
            if (spl[0] != "X0") continue;

            uint_t tid;
            std::string sid = spl[4];
            if (!begins_with(sid, "ID=")) {
                error("ill formed line, expected '# ID=...'");
                return 1;
            }

            sid = erase_begin(sid, "ID=");
            if (!from_string(sid, tid)) {
                error("could not read ID from '", sid, "'");
                return 1;
            }

            c_id.push_back(tid);
        }

        // print(i, ": ", c_id);

        for (uint_t k : range(c_id)) {
            uint_t j = where_first(incat.id == c_id[k]);
            if (j == npos) {
                error("could not find source ID=", c_id[k], " in input catalog");
                return 1;
            }

            outcat.x[j] = c_x0[k]-1;
            outcat.y[j] = c_y0[k]-1;
            outcat.n[j] = c_n[k];
            outcat.re[j] = c_re[k];
            outcat.pa[j] = c_pa[k];
            outcat.ell[j] = c_ell[k];

            vec2f cmod;
            std::string ofile_psf = base_out_dir+"models_conv/"+to_string(c_id[k])+".fits";
            fits::read(ofile_psf, cmod);
            outcat.flux[j] = total(cmod);
        }
    }

    fits::write(fit_sub_image, gimg, ghdr);

    fits::write_table(base_out_dir+"fluxes.fits", ftable(
        outcat.id, outcat.x, outcat.y, outcat.n, outcat.re, outcat.pa,
        outcat.ell, outcat.flux
    ));

    return 0;
}
