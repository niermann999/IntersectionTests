 cmake_minimum_required(VERSION 3.11)

project(Hyperplane VERSION 0.0 LANGUAGES CXX)

find_package(Eigen3 REQUIRED)
#set(Boost_NO_BOOST_CMAKE ON) 
find_package(Boost 1.69 REQUIRED COMPONENTS timer program_options unit_test_framework random)
find_package(Vc 1.4.1 REQUIRED) 

set(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/include)
set(CMAKE_BENCH_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/externals/benchmark/include)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/lib)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN_INCLUDE_DIRS})

# Define an interface target for the Vc library
add_library(Vc INTERFACE)
target_include_directories(Vc INTERFACE ${Vc_INCLUDE_DIR})
target_link_libraries(Vc INTERFACE libVc.a -L${Vc_LIB_DIR})

# Define an interface target for the core code
add_library(Core INTERFACE)
target_include_directories(Core INTERFACE ${CMAKE_INSTALL_INCLUDEDIR} ${CMAKE_SOURCE_DIR}/include/)
target_compile_features(Core INTERFACE cxx_std_17)

# General compile options
add_compile_options(
       -O3 --std=c++17
       #-L/usr/lib64/ -l:libVc.a
       #-Wall -Wextra -Wshadow -fstrict-aliasing -Wno-deprecated-declarations -Wstrict-aliasing
)

# Extra vectorization options
add_compile_options(
       #-march=native -msse -msse2 -msse3 -msse4 -mavx -mavx2
       #-msse -msse2
)

add_subdirectory(src)
add_subdirectory(externals)
