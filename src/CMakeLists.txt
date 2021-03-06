 macro(add_unittest _name)
  # automatically prefix the target name
  set(_target "UnitTest${_name}")
  add_executable(${_target} ${ARGN})

  # define required BOOST_TEST_... macros here to ensure consistent names
  target_compile_definitions(
    ${_target}
    PRIVATE "-DBOOST_TEST_DYN_LINK" "-DBOOST_TEST_MODULE=${_target}"
  )

  target_include_directories(
    ${_target}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  )

  target_link_libraries(
    ${_target}
    PRIVATE
    Core
    Boost::unit_test_framework Boost::timer Boost::random
    benchmark::benchmark
    Vc
    #Eigen3::Eigen
  )

  # register as unittest executable
  add_test(NAME ${_name} COMMAND ${_target})
endmacro()

macro(add_benchmark _name)
  # automatically prefix the target name
  set(_target "Benchmark${_name}")
  add_executable(${_target} ${ARGN})

  target_compile_definitions(
    ${_target}
    PRIVATE "-DBENCHMARK_ENABLE_GTEST_TESTS=OFF"
  )

  add_compile_options(-pthread)

  target_include_directories(
    ${_target}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
  )

  target_link_libraries(
    ${_target}
    PRIVATE
    Core
    benchmark::benchmark
    benchmark::benchmark_main
    Vc
    #Eigen3::Eigen
  )

  # register as unittest executable
  add_test(NAME ${_name} COMMAND ${_target})
endmacro()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_unittest(IntersectionTest IntersectionTests.cpp)
add_benchmark(EigenIntersection EigenBenchmark.cpp)
add_benchmark(EigenIntersection_wres EigenBenchmark_wres.cpp)
add_benchmark(VcVertIntersection VcVertBenchmark.cpp)
add_benchmark(VcVertIntersection_wres VcVertBenchmark_wres.cpp)
add_benchmark(VcHybridIntersection VcHybridBenchmark.cpp)
add_benchmark(VcHybridIntersection_wres VcHybridBenchmark_wres.cpp)
add_benchmark(VcHorizIntersection VcHorizBenchmark.cpp)
add_benchmark(VcHorizIntersection_wres VcHorizBenchmark_wres.cpp)
