cmake_minimum_required(VERSION 3.16)
project(epinetcpp)

find_package(GDAL REQUIRED)
find_package(CGAL REQUIRED)

include_directories(${GDAL_INCLUDE_DIR})
include_directories(${CGAL_INCLUDE_DIR})

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_STANDARD 17)
add_executable(epinetcpp nanoflann.hpp main.cpp main.h EpiMap.cpp EpiMap.h Common.h)
target_include_directories(epinetcpp PUBLIC "/usr/local/Cellar/gdal/2.4.4_5/include")
target_link_libraries(epinetcpp GDAL::GDAL CGAL::CGAL)
