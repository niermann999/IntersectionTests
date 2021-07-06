/**
 * author: joana.niermann@cern.ch
 **/
#include "fixtures.hpp"
#include "intersectors.hpp"

#include <benchmark/benchmark.h>

namespace vec_intr {

namespace g_bench {

//----------------------------------------------------//
// Vc Hybrid                                          //
//----------------------------------------------------//

//----------------------------------------------------Define Tests

//
// Use a gather on a structured data set then vectorize horizontaly 
// and save the results in a container
//
BENCHMARK_DEFINE_F(HybridSetup, intersectVcHybrid_wres)(benchmark::State& state) {
  using scalar_v = typename vector_v::vec_type;
  using vector_t = typename vector_v::type;

  aligned::vector<intersection<scalar_v, vector_t> > results;
  if (state.thread_index == 0) results.reserve(planes_struct.normals.size());
  //std::array<intersection<scalar_v, vector_t>, 10> results;

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    // Prevent compiler from optimizing away the loop
    benchmark::DoNotOptimize(results.data());

    vc_intersect_hybrid<vector_v>(ray_struct, planes_struct, results);

    //TODO: threadsafety
    results.clear();
  }
}

//----------------------------------------------------Run Tests

BENCHMARK_REGISTER_F(HybridSetup, intersectVcHybrid_wres)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  //->RangeMultiplier(n_surf_mult)->Range(n_surf_min, n_surf_steps*surf_step)
  ->Name("VcHybrid_wres")
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



