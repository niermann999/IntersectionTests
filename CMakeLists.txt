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
       -Wall -Wextra -Wshadow -fstrict-aliasing -Wno-deprecated-declarations -Wstrict-aliasing
)

# Set compiler flags manually (use -march=native as orientation)
if(COMPILE_SSE)
set(COMPILE_FLAGS -mmmx -msse -msse2 -msse3 -mssse3 -mno-sse4a -msse4.2 -msse4.1 -mno-avx -mno-avx2 -mno-3dnow -mno-avx512f -mno-avx512er -mno-avx512cd -mno-avx512pf -mno-avx512dq -mno-avx512bw -mno-avx512vl -mno-avx512ifma -mno-avx512vbmi -mno-avx5124fmaps -mno-avx5124vnniw -mno-avx512vbmi2 -mno-avx512vnni -mno-avx512bitalg -mno-avx512vpopcntdq)
elseif(COMPILE_AVX)
set(COMPILE_FLAGS -mavx -mavx2 -mno-avx512f -mno-avx512er -mno-avx512cd -mno-avx512pf  -mno-avx512dq -mno-avx512bw -mno-avx512vl -mno-avx512ifma -mno-avx512vbmi -mno-avx5124fmaps -mno-avx5124vnniw -mno-avx512vbmi2 -mno-avx512vnni -mno-avx512bitalg -mno-avx512vpopcntdq)
else()
set(COMPILE_FLAGS ${Vc_ARCHITECTURE_FLAGS})
endif()

message(STATUS "Using additional compiler flags: " ${COMPILE_FLAGS})
add_compile_options(
       ${COMPILE_FLAGS}
)

add_subdirectory(src)
add_subdirectory(externals)
