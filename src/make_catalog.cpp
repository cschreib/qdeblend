#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    std::string indir = argv[1];
    indir = file::directorize(indir);
    std::string fdbdir = argv[2];
    fdbdir = file::directorize(fdbdir);
    std::string id = argv[3];

    vec1s files = file::list_files(indir, "*allsub.fits");
    inplace_sort(files);

    auto db = read_filter_db(fdbdir+"/db.dat");
    auto fdb = read_filter_map(fdbdir+"fast.map");

    vec1s bands = erase_end(files, "-allsub.fits");
    vec1f lambda = replicate(fnan, bands.size());
    vec1f lambda_low = replicate(fnan, bands.size());
    vec1f lambda_up = replicate(fnan, bands.size());
    vec1u fbands = replicate(npos, bands.size());
    vec1f zp = replicate(1.0, bands.size());

    if (!file::exists("zero-points.txt")) {
        ascii::write_table("zero-points.txt", bands, zp);
    } else {
        vec1s tbands;
        vec1f tzp;
        ascii::read_table("zero-points.txt", tbands, tzp);
        vec1u id1, id2;
        match(bands, tbands, id1, id2);
        zp[id1] = tzp[id2];
    }

    for (uint_t b : range(bands)) {
        auto& s = bands[b];

        filter_t fil;
        if (get_filter(db, s, fil)) {
            lambda[b] = fil.rlam;
            lambda_low[b] = fil.lam[where_first(fil.res > 0.5*max(fil.res))];
            lambda_up[b] = fil.lam[where_last(fil.res > 0.5*max(fil.res))];
        }

        uint_t ifast;
        if (get_filter_id(fdb, s, ifast)) {
            fbands[b] = ifast;
        } else {
            warning("could not find FAST ID for band ", bands[b]);
        }
    }

    vec1s sid;
    vec1d ra, dec;
    vec2f flux, flux_err;

    for (uint_t b : range(bands)) {
        vec1s model;
        vec1d tflx, terr, terr_sim;
        fits::read_table(indir+files[b], "flux", tflx, "flux_err", terr, "flux_err_sim", terr_sim,
            "model", model);

        if (count(is_finite(terr_sim)) == 0) {
            warning("band '", bands[b], "' is missing simulations");
        }

        terr = max(terr, terr_sim);

        // Apply zero point correction
        tflx /= zp[b];
        terr /= zp[b];

        // Clamp error to 10% of flux (max S/N=10)
        terr = max(terr, 0.1*tflx);

        vec1s tsid = file::remove_extension(file::get_basename(model));

        if (sid.empty()) {
            sid = tsid;

            flux = replicate(fnan, sid.size(), bands.size());
            flux_err = replicate(fnan, sid.size(), bands.size());

            ra = replicate(dnan, sid.size());
            dec = replicate(dnan, sid.size());

            for (uint_t i : range(sid)) {
                std::string ff = indir+"imfit/models/"+sid[i]+".fits";

                fits::input_image iimg(ff);
                vec2d data;
                iimg.read(data);
                vec1u mid = mult_ids(data, max_id(data));

                astro::xy2ad(astro::wcs(iimg.read_header()), mid[1]+1, mid[0]+1, ra[i], dec[i]);
            }

            flux(_,b) = tflx;
            flux_err(_,b) = terr;
        } else {
            vec1u id1, id2;
            match(sid, tsid, id1, id2);
            id2 = complement(tsid, id2);

            append(sid, tsid[id2]);
            append(ra, replicate(dnan, id2.size()));
            append(dec, replicate(dnan, id2.size()));
            append<0>(flux, replicate(fnan, id2.size(), flux.dims[1]));
            append<0>(flux_err, replicate(fnan, id2.size(), flux.dims[1]));

            match(sid, tsid, id1, id2);
            flux(id1,b) = tflx[id2];
            flux_err(id1,b) = terr[id2];
        }
    }

    // In FITS format
    fits::write_table("catalog.fits", ftable(
        sid, ra, dec, flux, flux_err, bands, lambda, lambda_low, lambda_up
    ));

    // In FAST-like format with all bands
    vec1s hdr_fast = {"id", "z_spec", "ra", "dec"};
    for (uint_t b : range(bands)) {
        if (fbands[b] != npos) {
            hdr_fast.push_back("F"+to_string(fbands[b]));
            hdr_fast.push_back("E"+to_string(fbands[b]));
        }
    }

    // In true FAST format with only the bands FAST knows
    vec1u idf = where(fbands != npos);
    vec1u ids = where(sid == id);
    ascii::output_format opts;
    opts.header = hdr_fast;
    ascii::write_table("fast.cat", opts,
        vec1u{1}, vec1d{-1.0}, to_string_vector(format::precision(ra[ids], 10)),
        to_string_vector(format::precision(dec[ids], 10)),
        ascii::columns(idf.size(), flux(ids,idf), flux_err(ids,idf))
    );

    ascii::write_table("fast_all.cat", opts,
        sid, replicate(-1.0, sid.size()), to_string_vector(format::precision(ra, 10)),
        to_string_vector(format::precision(dec, 10)),
        ascii::columns(idf.size(), flux(_,idf), flux_err(_,idf))
    );

    return 0;
}
