# IntersectionTests
Toy example of planar surface intersection test to optimize vectorization

Depends on:
- Eigen3 
- Boost 1.69 (timer program_options unit_test_framework)
- Vc 1.4.1

To build the project do:
git clone https://github.com/niermann999/IntersectionTests.git
mkdir IntesectionTests/build
cd IntersectionTests/build
git submodule update --init --recursive

Run cmake, e.g.
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER="g++" -DCMAKE_CXX_FLAGS="-O3 -march=native -msse -msse2 -msse3 -msse4 -mmmx -m3dnow -mavx -mavx2 -mfma" ..
make -j4

Some caution is still needed with compiler flags when using cmake
(-march=native -msse -msse2 -msse3 -msse4 -mavx -mavx2)
