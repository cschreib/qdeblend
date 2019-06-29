#include <vif.hpp>
#include <vif/astro/ds9.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    std::string bad_file = argv[2];

    fits::header hdr;
    astro::wcs w; {
        fits::input_image iimg(argv[1]);
        hdr = iimg.read_header();
        w = astro::wcs(hdr);
    }

    // Mask bad pixels
    vec2b bad(w.dims[0], w.dims[1]);
    ds9::mask_regions(bad_file, w, bad);

    fits::output_image oimg(argv[3]);
    oimg.write(bad);
    oimg.write_header(hdr);

    return 0;
}
