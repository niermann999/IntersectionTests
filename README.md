# IntersectionTests
Toy example of planar surface intersection test to optimize vectorization

Depends on:
- Eigen3 
- Boost 1.69 (timer program_options unit_test_framework)
- Vc 1.4.1

For now you need to install google benchmark at "IntersectionTests/benchmark" by
hand, as described in their doc.

Some caution is still needed with compiler flags when using cmake
(-march=native -msse -msse2 -msse3 -msse4 -mavx -mavx2)
