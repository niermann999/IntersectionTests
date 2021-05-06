/**
 * author: joana.niermann@cern.ch
 **/
#include "fixtures.hpp"
#include "intersectors.hpp"

#include <benchmark/benchmark.h>

namespace vec_intr {

namespace g_bench {
  
//----------------------------------------------------//
// Vc Vertical                                        //
//----------------------------------------------------//

//----------------------------------------------------Define Tests
//
// Use simdArray on unchanged data set
//
BENCHMARK_DEFINE_F(VertSetup, intersectVcVert)(benchmark::State& state) {

  for (auto _: state) {
    for (auto &plane : planes) {
      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        vc_intersect_vert<vector_s>(ray, plane)
      );
      // Prevent compiler from optimizing away the loop
      benchmark::ClobberMemory();
    }
  }
}

//
// Use simdArray on unchanged data set and save the results in a container
//
BENCHMARK_DEFINE_F(VertSetup, intersectVcVert_wres)(benchmark::State& state) {
  using scalar_t = typename vector_s::scalar_type;

  aligned::vector<intersection<scalar_t, Vc::SimdArray<scalar_t, 4> > > results;
  if (state.thread_index == 0) results.reserve(planes.size());

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    // Allow vector data to be clobbered in memory
    benchmark::DoNotOptimize(results.data());

    vc_intersect_vert<vector_s>(ray, planes, results);

    // Prevent compiler from optimizing away the loop
    benchmark::ClobberMemory();
    //TODO: Check threadsafety
    results.clear();
  }
}

//----------------------------------------------------Run Tests

BENCHMARK_REGISTER_F(VertSetup, intersectVcVert)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("VcVert")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef NO_MULTI_THREAD
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(VertSetup, intersectVcVert_wres)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("VcVert_wres")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef NO_MULTI_THREAD
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif


BENCHMARK_MAIN();

} // namespace g_bench

} //namespace vec_intr



