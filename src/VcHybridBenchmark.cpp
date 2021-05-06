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
//
BENCHMARK_DEFINE_F(HybridSetup, intersectVcHybrid)(benchmark::State& state) {
  using scalar_t = typename vector_v::scalar_type;
  using scalar_v = typename vector_v::vec_type;

  for (auto _: state) {
    for (Index_v i(scalar_v::IndexesFromZero()); (i < Index_v(planes_struct.points.size())).isFull(); i += Index_v(scalar_v::Size)) {
      vector_v pl_normal_strc, pl_point_strc;

      scalar_v pns_x = planes_struct.normals[i][&Vector3<scalar_t>::x];
      scalar_v pns_y = planes_struct.normals[i][&Vector3<scalar_t>::y];
      scalar_v pns_z = planes_struct.normals[i][&Vector3<scalar_t>::z];
      pl_normal_strc.obj = {.x = pns_x, .y = pns_y, .z = pns_z};

      scalar_v pps_x = planes_struct.points[i][&Vector3<scalar_t>::x];
      scalar_v pps_y = planes_struct.points[i][&Vector3<scalar_t>::y];
      scalar_v pps_z = planes_struct.points[i][&Vector3<scalar_t>::z];
      pl_point_strc.obj = {.x = pps_x, .y = pps_y, .z = pps_z};
      plane_data<vector_v> planes_strcts = {.normals = pl_normal_strc, .points = pl_point_strc};

      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        vc_intersect_hybrid<vector_v>(ray_struct, planes_strcts)
      );
      // Prevent compiler from optimizing away the loop
      benchmark::ClobberMemory();
    }
  }
}

//
// Use a gather on a structured data set then vectorize horizontaly 
// and save the results in a container
//
BENCHMARK_DEFINE_F(HybridSetup, intersectVcHybrid_wres)(benchmark::State& state) {
  using scalar_v = typename vector_v::vec_type;
  using vector_t = typename vector_v::type;

  aligned::vector<intersection<scalar_v, vector_t> > results;
  if (state.thread_index == 0) results.reserve(planes_struct.normals.size());

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    // Allow vector data to be clobbered in memory
    benchmark::DoNotOptimize(results.data());

    vc_intersect_hybrid<vector_v>(ray_struct, planes_struct, results);

    // Prevent compiler from optimizing away the loop
    benchmark::ClobberMemory();
    //TODO: threadsafety
    results.clear();
  }
}

//----------------------------------------------------Run Tests

BENCHMARK_REGISTER_F(HybridSetup, intersectVcHybrid)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("VcHybrid")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef NO_MULTI_THREAD
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(HybridSetup, intersectVcHybrid_wres)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("VcHybrid_wres")
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



