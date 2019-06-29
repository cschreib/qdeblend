#include <vif.hpp>
#include <vif/astro/ds9.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    std::string bad_file = argv[2];

    fits::image iimg(argv[1]);
    fits::header hdr = iimg.read_header();
    astro::wcs w(hdr);

    vec2b bad;
    iimg.read(bad);
    ds9::mask_regions(bad_file, w, bad);

    iimg.update(bad);

    return 0;
}
