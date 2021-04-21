 macro(add_unittest _name)
  # automatically prefix the target name
  set(_target "UnitTest${_name}")
  add_executable(${_target} ${ARGN})

  target_compile_features(${_target} PUBLIC cxx_std_17)

  # define required BOOST_TEST_... macros here to ensure consistent names
  target_compile_definitions(
    ${_target}
    PRIVATE "-DBOOST_TEST_DYN_LINK" "-DBOOST_TEST_MODULE=${_target}")
 
  target_include_directories(
   ${_target}
   SYSTEM PUBLIC ${EIGEN_INCLUDE_DIRS}
  )

  target_include_directories(${_target} 
   PRIVATE ${CMAKE_INSTALL_INCLUDEDIR}
  )

  target_include_directories(
    ${_target}
    PRIVATE Core ${CMAKE_CURRENT_SOURCE_DIR})

  target_link_libraries(
    ${_target}
    PRIVATE
      Boost::unit_test_framework Boost::timer
      benchmark::benchmark
      libVc.a)
  # register as unittest executable
  add_test(NAME ${_name} COMMAND ${_target})
endmacro()

macro(add_benchmark _name)
  # automatically prefix the target name
  set(_target "Benchmark${_name}")
  add_executable(${_target} ${ARGN})

  target_compile_features(${_target} PUBLIC cxx_std_17)

  target_include_directories(${_target} 
   PRIVATE ${CMAKE_BENCH_INCLUDEDIR}
  )
 
  target_include_directories(
   ${_target}
   SYSTEM PUBLIC ${EIGEN_INCLUDE_DIRS}
  )

  target_include_directories(${_target} 
   PRIVATE ${CMAKE_INSTALL_INCLUDEDIR}
  )

  target_include_directories(
    ${_target}
    PRIVATE Core ${CMAKE_CURRENT_SOURCE_DIR})

  target_link_libraries(
    ${_target}
    PRIVATE
      benchmark::benchmark
      benchmark::benchmark_main
      libVc.a)
  # register as unittest executable
  add_test(NAME ${_name} COMMAND ${_target})
endmacro()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

add_unittest(IntersectionTest IntersectionTests.cpp)
add_benchmark(IntersectionBench IntersectionBenchmark.cpp)
