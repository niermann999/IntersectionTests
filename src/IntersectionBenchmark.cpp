 /**
 * author: asalzburger@gmail.com
 **/
#include <types.hpp>
#include <intersectors.hpp>

#include <array>
#include <ctime>            // std::time
#include <chrono>
#include <iostream>
#include <limits>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

#include <benchmark/benchmark.h>

namespace vec_intr {

namespace g_bench {

static void BM_StringCreation(benchmark::State& state) {
  for (auto _ : state)
    std::string empty_string;
}
// Register the function as a benchmark
BENCHMARK(BM_StringCreation);

// Define another benchmark
static void BM_StringCopy(benchmark::State& state) {
  std::string x = "hello";
  for (auto _ : state)
    std::string copy(x);
}
BENCHMARK(BM_StringCopy);

BENCHMARK_MAIN();


} // namespace g_bench

} //namespace vec_intr


