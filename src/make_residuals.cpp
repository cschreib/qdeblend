#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    vec1s files = file::list_files("../", "*-allsub-sub.fits");
    files = replace(files, "-allsub-sub.fits", "");

    file::mkdir("residuals/before");
    file::mkdir("residuals/after");

    std::string id = argv[1];

    for (auto& b : files) {
        std::string filename = "../"+b+"-allsub-sub.fits";
        fits::input_image iimg(filename);

        bool found = false;
        for (uint_t i : range(1, iimg.hdu_count())) {
            iimg.reach_hdu(i);
            std::string extname;
            if (iimg.read_keyword("EXTNAME", extname) && extname == "RES_"+id) {
                found = true;
                break;
            }
        }

        if (!found) {
            warning("could not find source '", id, "' in ", filename);
            continue;
        }

        vec2d img;
        iimg.read(img);
        fits::header hdr = iimg.read_header();

        fits::write("residuals/before/"+b+".fits", vec2f{img}, hdr);
        file::copy("../"+b+".fits", "residuals/after/"+b+".fits");
    }

    return 0;
}
