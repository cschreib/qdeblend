#include <vif.hpp>
#include <vif/astro/ds9.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    vec1s files = {"acs-f435w", "acs-f606w", "acs-f775w", "acs-f814w", "acs-f850lp",
        "wfc3-f105w", "wfc3-f125w", "wfc3-f140w", "wfc3-f160w"};
    files += ".fits";

    std::string flux_image = "wfc3-f160w.fits";

    double deblend_threshold = 0.4;
    double detect_threshold = 1.0;
    double conv_fwhm = 0.2;
    double min_snr = 10.0;
    uint_t min_area = 10;
    read_args(argc, argv, arg_list(files, flux_image,
        deblend_threshold, detect_threshold, conv_fwhm, min_area, min_snr
    ));

    vec3f imgs;
    vec1f rms;
    double stack_rms = 0;
    vec2f stack;
    vec2f wstack;
    fits::header ghdr;
    double aspix = 0.06;

    for (std::string& f : files) {
        if (!file::exists(f)) continue;

        vec2d img;
        fits::header hdr;
        fits::read(f, img, hdr);
        astro::wcs w(hdr);

        if (ghdr.empty()) ghdr = hdr;

        // From DS9 region files
        std::string base = file::remove_extension(f);
        std::string bad_file = base+"-bad.reg";
        if (!file::exists(bad_file)) {
            bad_file = base+"_bad.reg";
        }
        if (!file::exists(bad_file)) {
            bad_file = "";
        }

        if (!bad_file.empty()) {
            vec2b bad_mask(img.dims);
            ds9::mask_regions(bad_file, w, bad_mask);
            img[where(bad_mask)] = dnan;
        }

        if (fraction_of(is_finite(img) || abs(img) > 0) < 0.6) continue;

        rms.push_back(1.48*mad(img));

        vec2d wei = replicate(1.0/sqr(rms.back()), img.dims);
        vec1u idb = where(!is_finite(img));
        img[idb] = 0; wei[idb] = 0;

        if (stack.empty()) {
            stack = img*wei;
            wstack = wei;
        } else {
            stack += img*wei;
            wstack += wei;
        }

        stack_rms += sqr(rms.back()/sqr(rms.back()));

        append<0>(imgs, reform(img/rms.back(), 1, img.dims));
    }

    astro::wcs gw(ghdr);

    double tw = total(1/sqr(rms));
    stack_rms = sqrt(stack_rms)/tw;
    stack /= wstack;
    fits::write("det_stack.fits", stack, ghdr);

    append<0>(imgs, reform(stack/stack_rms, 1, stack.dims));

    // Build detection image
    vec2f det = partial_max(0, imgs);

    double conv = conv_fwhm/aspix/2.335; {
        vec1u idb = where(!is_finite(det));
        det[idb] = 0;
        uint_t hpix = 5*conv;
        vec2d kernel = gaussian_profile({{2*hpix+1, 2*hpix+1}}, conv, hpix, hpix);
        det = convolve2d(det, kernel);
        det[idb] = fnan;
    }

    fits::write("det_snr.fits", det, ghdr);

    // Perform segmentation
    segment_deblend_params sdp;
    sdp.detect_threshold = detect_threshold + median(det);
    sdp.deblend_threshold = deblend_threshold/aspix;
    sdp.min_area = min_area;

    segment_deblend_output segments;
    vec2u seg = segment_deblend(det, segments, sdp);

    // Read fluxes on image
    vec2d himg;
    fits::read(flux_image, himg);

    vec1u idbg = where(seg == 0);
    himg -= median(himg[idbg]);

    double hrms = 1.48*mad(himg[idbg]);

    vec1f hflx(segments.id.size());
    vec1f hflx_err(segments.id.size());
    foreach_segment(seg, segments.origin, [&](uint_t s, vec1u ids) {
        hflx[s] = total(himg[ids]);
        hflx_err[s] = sqrt(ids.size())*hrms;
    });

    // Remove too low S/N objects.
    vec1u idls = where(hflx/hflx_err < min_snr);
    foreach_segment(seg, segments.origin[idls], [&](uint_t s, vec1u ids) {
        seg[ids] = 0;
    });

    inplace_remove(segments.id, idls);
    inplace_remove(segments.area, idls);
    inplace_remove(segments.px, idls);
    inplace_remove(segments.py, idls);
    inplace_remove(hflx, idls);
    inplace_remove(hflx_err, idls);

    // Merge segments (if any).
    std::string manual_file = "merged_segments.reg";
    if (file::exists(manual_file)) {
        vec<1,ds9::region> reg;
        ds9::read_regions_physical(manual_file, gw, reg);

        for (uint_t i : range(reg)) {
            vec2b mask(seg.dims);
            ds9::mask_region(reg[i], mask);

            vec1u idm = where(mask);
            vec1u idu = unique_values(seg[idm]);
            idu = idu[where(idu != 0)];

            if (idu.size() > 1) {
                print("merged ", idu.size(), " segments");
                for (uint_t u : range(1, idu.size())) {
                    vec1u ids = where(seg == idu[u]);
                    seg[ids] = idu[0];
                    idu[u] = where_first(segments.id == idu[u]);
                }

                uint_t im = where_first(segments.id == idu[0]);
                idu[0] = im;

                segments.px[im] = total(segments.px[idu]*segments.area[idu])/total(segments.area[idu]);
                segments.py[im] = total(segments.py[idu]*segments.area[idu])/total(segments.area[idu]);
                segments.area[im] = total(segments.area[idu]);
                hflx[im] = total(hflx[idu]);
                hflx_err[im] = sqrt(total(sqr(hflx_err[idu])));

                idu = idu[1-_];
                inplace_remove(segments.id, idu);
                inplace_remove(segments.px, idu);
                inplace_remove(segments.py, idu);
                inplace_remove(segments.area, idu);
                inplace_remove(hflx, idu);
                inplace_remove(hflx_err, idu);
            }
        }
    }


    // Add manual detections (if any).
    manual_file = "det_manual.reg";
    if (file::exists(manual_file)) {
        vec<1,ds9::region> reg;
        ds9::read_regions_physical(manual_file, gw, reg);

        vec1s sid(reg.size());
        for (uint_t i : range(reg)) sid[i] = reg[i].text;

        vec1s uid = unique_values(sid);
        for (std::string s : uid) {
            vec1u idl = where(sid == s);
            vec2b mask(seg.dims);
            for (uint_t i : idl) {
                ds9::mask_region(reg[i], mask);
            }

            vec1u ids = where(mask);
            uint_t seg_id = max(segments.id)+1;
            seg[ids] = seg_id;

            segments.id.push_back(seg_id);
            segments.px.push_back(reg[idl[0]].params[0]);
            segments.py.push_back(reg[idl[0]].params[1]);
            segments.area.push_back(ids.size());
            hflx.push_back(total(himg[ids]));
            hflx_err.push_back(sqrt(ids.size())*hrms);
        }
    }

    fits::write("det_seg.fits", seg, ghdr);

    vec1d ra, dec;
    astro::xy2ad(gw, segments.px+1, segments.py+1, ra, dec);

    fits::write_table("det_cat.fits",
        "id", segments.id, "area", segments.area, "x", segments.px, "y", segments.py,
        "ra", ra, "dec", dec, "flux", hflx, "flux_err", hflx_err
    );

    return 0;
}
