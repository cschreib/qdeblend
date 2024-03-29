cmake_minimum_required(VERSION 2.6)
project(speclib-superbuild C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected (-DCMAKE_BUILD_TYPE=...)")
    message(STATUS "Default to Release.")
    set(CMAKE_BUILD_TYPE "Release")
endif()

include(ExternalProject)

# Set default install dir to "here"
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "No install directory selected (-DCMAKE_INSTALL_PREFIX=...)")
    message(STATUS "Default to ${PROJECT_SOURCE_DIR}/bin")
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "default install path" FORCE)
endif()

# Configure and build cfitsio
ExternalProject_Add(cfitsio
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/cfitsio
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/cfitsio
    GIT_REPOSITORY https://github.com/HEASARC/cfitsio.git
    GIT_TAG ff6bdb0c2f891ae6709f6bb2378e24e005b7c92b
    CMAKE_ARGS
        -DBUILD_SHARED_LIBS=OFF
        -DUSE_PTHREADS=ON
        -DBUILD_TOOLS=OFF
        -DBUILD_TESTS=OFF
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build WCSlib
ExternalProject_Add(wcslib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/wcslib
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/wcslib
    URL ftp://ftp.atnf.csiro.au/pub/software/wcslib/wcslib.tar.bz2
    PATCH_COMMAND patch ${PROJECT_SOURCE_DIR}/external/wcslib/makedefs.in ${PROJECT_SOURCE_DIR}/patches/wcslib_shrlib.patch
    CONFIGURE_COMMAND ${PROJECT_SOURCE_DIR}/external/wcslib/configure
        --prefix=${PROJECT_SOURCE_DIR}/external
        --without-cfitsio --without-pgplot --disable-fortran --disable-utils
        CC=${CMAKE_C_COMPILER}
    BUILD_COMMAND make install
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1)

# Configure and build vif
ExternalProject_Add(vif
    DEPENDS cfitsio wcslib
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/external/vif
    DOWNLOAD_DIR ${PROJECT_SOURCE_DIR}/external/vif
    GIT_REPOSITORY https://github.com/cschreib/vif.git
    GIT_TAG 61653ddd286f9bc49ceb594aeb089dc48c71d272
    CMAKE_ARGS
        -DNO_REFLECTION=ON
        -DNO_FFTW=OFF
        -DNO_LAPACK=ON
        -DNO_GSL=ON
        -DVIF_INPLACE_BUILD=ON
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DWCSLIB_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=Release)

# Configure and build sub-projects
ExternalProject_Add(qdeblend
    DEPENDS vif
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/src
    CMAKE_ARGS
        -DCFITSIO_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DWCSLIB_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DVIF_ROOT_DIR=${PROJECT_SOURCE_DIR}/external
        -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}
        -DCMAKE_MODULE_PATH=${PROJECT_SOURCE_DIR}/external/vif/cmake
        -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
        -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
