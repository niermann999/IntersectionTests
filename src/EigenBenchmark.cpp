/**
 * author: joana.niermann@cern.ch
 **/
#include "fixtures.hpp"
#include "intersectors.hpp"

#include <benchmark/benchmark.h>

namespace vec_intr {

namespace g_bench {


//----------------------------------------------------//
// Eigen                                              //
//----------------------------------------------------//

//----------------------------------------------------Define Tests

BENCHMARK_DEFINE_F(VertSetup, intersectEigen4D)(benchmark::State& state) {

  for (auto _: state) {
    for (auto &plane : planes) {
      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        eig_intersect_4D<vector_s>(ray, plane)
      );
      // Prevent compiler from optimizing away the loop
      benchmark::ClobberMemory();
    }
  }
}

BENCHMARK_DEFINE_F(VertSetup, intersectEigen4D_wres)(benchmark::State& state) {
  using scalar_t = typename vector_s::scalar_type;
  using vector_t = typename vector_s::type;

  aligned::vector<intersection<scalar_t, vector_t> > results;
  if (state.thread_index == 0) results.reserve(planes.size());

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    // Allow vector data to be clobbered in memory
    benchmark::DoNotOptimize(results.data());

    eig_intersect_4D<vector_s>(ray, planes, results);

    // Prevent compiler from optimizing away the loop
    benchmark::ClobberMemory();
    //TODO: Check threadsafety
    results.clear();
  }
}

//----------------------------------------------------Run Tests

BENCHMARK_REGISTER_F(VertSetup, intersectEigen4D)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("Eigen4D")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef MULTI_THREAD
  ->Threads(MULTI_THREAD);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(VertSetup, intersectEigen4D_wres)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("Eigen4D_wres")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef MULTI_THREAD
  ->Threads(MULTI_THREAD);
  #else
  ->ThreadPerCpu();
  #endif


BENCHMARK_MAIN();

} // namespace g_bench

} //namespace vec_intr



