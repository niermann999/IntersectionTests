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
      // Prevent compiler from optimizing away the loop
      benchmark::DoNotOptimize(
        vc_intersect_vert<vector_s>(ray, plane)
      );
    }
  }
}

//----------------------------------------------------Run Tests

BENCHMARK_REGISTER_F(VertSetup, intersectVcVert)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  //->RangeMultiplier(n_surf_mult)->Range(n_surf_min, n_surf_steps*surf_step)
  ->Name("VcVert")
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



