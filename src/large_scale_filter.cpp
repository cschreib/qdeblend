#include <vif.hpp>

using namespace vif;
using namespace vif::astro;

int vif_main(int argc, char* argv[]) {
    std::string image;
    std::string mask;
    double size = 0;

    read_args(argc, argv, arg_list(image, mask, size));

    vec2d img;
    fits::header hdr;
    fits::read(image, img, hdr);
    vec2b bad;
    fits::read(mask, bad);

    double aspix = 0;
    get_pixel_size(astro::wcs(hdr), aspix);

    uint_t isize = round(size/aspix)/2;
    if (2*isize+1 > img.dims[0] || 2*isize+1 > img.dims[1]) {
        warning("box is larger than image, not doing anything");
        return 0;
    }

    vec2d nimg = img;
    img[where(bad)] = dnan;
    for (uint_t iy : range(img.dims[0]))
    for (uint_t ix : range(img.dims[1])) {
        uint_t y0 = iy > isize ?               iy - isize : 0;
        uint_t y1 = iy < img.dims[0] - isize ? iy + isize : img.dims[0] - 1;
        uint_t x0 = ix > isize ?               ix - isize : 0;
        uint_t x1 = ix < img.dims[1] - isize ? ix + isize : img.dims[1] - 1;

        nimg.safe(iy,ix) -= median(img(y0-_-y1,x0-_-x1));
    }

    fits::write(image, nimg, hdr);


    return 0;
}
