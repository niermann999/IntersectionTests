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

      // Prevent compiler from optimizing away the loop
      benchmark::DoNotOptimize(
        vc_intersect_horiz<vector_v>(ray_struct, planes_strcts)
      );
    }
  }
}

//----------------------------------------------------Run Tests

BENCHMARK_REGISTER_F(HybridSetup, intersectVcHybrid)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  //->RangeMultiplier(n_surf_mult)->Range(n_surf_min, n_surf_steps*surf_step)
  ->Name("VcHybrid")
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



