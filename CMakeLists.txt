 cmake_minimum_required(VERSION 3.11)

project(Hyperplane VERSION 0.0 LANGUAGES CXX)

find_package(Eigen3 REQUIRED)
set(Boost_NO_BOOST_CMAKE ON) 
find_package(Boost 1.69 REQUIRED COMPONENTS timer program_options unit_test_framework)

find_package(Vc 1.4.1 REQUIRED) 

set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

add_subdirectory(src)
