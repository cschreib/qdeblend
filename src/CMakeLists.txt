cmake_minimum_required(VERSION 2.6)
project(xshooter C CXX)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(NO_REFLECTION ON)
set(NO_FFTW ON)
set(NO_GSL ON)
set(NO_LAPACK ON)

find_package(vif REQUIRED)

# Include dependencies
message(STATUS ${VIF_INCLUDE_DIRS})
include_directories(${VIF_INCLUDE_DIRS})

# Build tools
add_executable(add_mask add_mask.cpp)
target_link_libraries(add_mask ${VIF_LIBRARIES})
install(TARGETS add_mask DESTINATION bin)

add_executable(go_imfit go_imfit.cpp)
target_link_libraries(go_imfit ${VIF_LIBRARIES})
install(TARGETS go_imfit DESTINATION bin)

add_executable(lage_scale_filter lage_scale_filter.cpp)
target_link_libraries(lage_scale_filter ${VIF_LIBRARIES})
install(TARGETS lage_scale_filter DESTINATION bin)

add_executable(make_detect make_detect.cpp)
target_link_libraries(make_detect ${VIF_LIBRARIES})
install(TARGETS make_detect DESTINATION bin)

add_executable(make_mask make_mask.cpp)
target_link_libraries(make_mask ${VIF_LIBRARIES})
install(TARGETS make_mask DESTINATION bin)

add_executable(multi_subtract multi_subtract.cpp)
target_link_libraries(multi_subtract ${VIF_LIBRARIES})
install(TARGETS multi_subtract DESTINATION bin)
