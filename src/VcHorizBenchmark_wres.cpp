/**
 * author: joana.niermann@cern.ch
 **/
#include "fixtures.hpp"
#include "intersectors.hpp"

#include <benchmark/benchmark.h>
#ifdef DEBUG
#include <chrono>
#include <ctime>
#endif

namespace vec_intr {

namespace g_bench {


#ifdef DEBUG
// Measure execution time
using clock = std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using unit_ms = std::chrono::milliseconds;
#endif

//----------------------------------------------------//
// Vc Horizontal                                      //
//----------------------------------------------------//

//----------------------------------------------------Define Tests

//
// Use a horizontal vectorization and data set, save the results in a container
//
BENCHMARK_DEFINE_F(HorizSetup, intersectVcHoriz_wres)(benchmark::State& state) {
  using scalar_v = typename vector_v::vec_type;
  using vector_t = typename vector_v::type;

  aligned::vector<intersection<scalar_v, vector_t> > results;
  if (state.thread_index == 0) results.reserve(planes_hor.size());

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    // Prevent compiler from optimizing away the loop
    benchmark::DoNotOptimize(results.data());

    vc_intersect_horiz<vector_v>(ray_hor, planes_hor, results);

    //TODO: threadsafety
    results.clear();
  }
}

//----------------------------------------------------Run Tests

BENCHMARK_REGISTER_F(HorizSetup, intersectVcHoriz_wres)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  //->RangeMultiplier(n_surf_mult)->Range(n_surf_min, n_surf_steps*surf_step)
  ->Name("VcHoriz_wres")
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



