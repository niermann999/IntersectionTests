 /**
 * author: joana.niermann@cern.ch
 * TODO: isolate benchmarks from each other
 **/
#include "types.hpp"
#include "intersectors.hpp"

#ifdef DEBUG
#include <chrono>
#include <ctime>
#endif
#include <iostream>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

#include <benchmark/benchmark.h>

namespace vec_intr {

namespace g_bench {

// Iterations of a benchmark
//constexpr size_t gbench_test_itrs = 10000;
// Repetitions of a benchmark
constexpr size_t gbench_test_repts = 10;
// Number of rand. gen. surfaces to intersect
constexpr size_t surf_step    = 5;
constexpr size_t n_surf_steps = 100;
#ifdef NO_MULTI_THREAD
constexpr size_t nThreads  = 1;
#endif

// Make sure the memory layout is compatible with Vc Vectors and set corresponding LA wrappers as types
static_assert(data_trait<Vector4_s>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");
static_assert(data_trait<Vector3_v>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");

using vector_s = data_trait<Vector4_s>::type;
using vector_v = data_trait<Vector3_v>::type; // Contains a vector in every coordinate

// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
using generator_t = boost::minstd_rand;
// Define a uniform random number distribution which produces "double"
// values between 0 and 1 (0 inclusive, 1 exclusive).
using rand_t = boost::variate_generator<generator_t&, boost::uniform_real<Scalar> >;

#ifdef DEBUG
// Measure execution time
using clock = std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using unit_ms = std::chrono::milliseconds;
#endif

//----------------------------------------------------Fill Data Types

generator_t generator(42);

boost::uniform_real<Scalar> uni_dist(0,1000);
rand_t uni(generator, uni_dist);

// Vertical data as in detray
class VertSetup : public benchmark::Fixture {

    public:

    vector_s ray_dir, ray_point;
    vector_s pl_normals, pl_points;

    ray_data<vector_s> ray;
    aligned::vector<plane_data<vector_s> > planes;

    #ifdef DEBUG
    // manual time measurement using chrono
    unit_ms duration;
    #endif

    ~VertSetup() = default;
    void SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      ray_dir.obj   = vector_s::type::Random();
      ray_point.obj = vector_s::type::Random();

      planes.reserve(state.range(0));
      for (size_t i = 0; i < state.range(0); i++) {
        pl_normals = {.obj = vector_s::type::Random()};
        pl_points  = {.obj = vector_s::type::Random()};

        planes.push_back({.normals = pl_normals, .points = pl_points});
      }

      // vertical data containers
      ray = {.direction = std::move(ray_dir), 
             .point     = std::move(ray_point)};
    }

    void TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes.clear();
      }
    }
    
};

// Vertical data contained in structs
class HybridSetup : public benchmark::Fixture {

    public:

    vector_v ray_dir_hor, ray_point_hor;
    aligned::vector<Vector3<Scalar> > pl_points_struct, pl_normals_struct;

    ray_data<vector_v> ray_struct;
    plane_data<aligned::vector<Vector3<Scalar> > > planes_struct;

    ray_data<vector_v> ray_hor;
    plane_data<aligned::vector<vector_v> > planes_hor;

    #ifdef DEBUG
    // manual time measurement using chrono
    unit_ms duration;
    #endif

    ~HybridSetup() = default;
    void SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      // AoS data
      pl_normals_struct.reserve(state.range(0));
      pl_points_struct.reserve(state.range(0));
      for (size_t i = 0; i < state.range(0); i++) {
        pl_normals_struct.push_back({.x=uni(), .y=uni(), .z=uni()});
        pl_points_struct.push_back( {.x=uni(), .y=uni(), .z=uni()});
      }

      // Horizontal data (interleaved)
      ray_dir_hor.obj   = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
      ray_point_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};

      // AoS container
      ray_struct    = {.direction = std::move(ray_dir_hor), 
                       .point     = std::move(ray_point_hor)}; 
      planes_struct = {.normals = std::move(pl_normals_struct), 
                       .points  = std::move(pl_points_struct)};
    }

    void TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes_struct.normals.clear();
        planes_struct.points.clear();
      }
    }
    
};

// Horizontal data as interleaved horizonral vectors
class HorizSetup : public benchmark::Fixture {

    public:

    vector_v ray_dir_hor, ray_point_hor;
    vector_v pl_normals_hor, pl_points_hor;

    ray_data<vector_v> ray_hor;
    aligned::vector<plane_data<vector_v> > planes_hor;

    #ifdef DEBUG
    // manual time measurement using chrono
    unit_ms duration;
    #endif

    ~HorizSetup() = default;
    void SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      // Horizontal data (interleaved)
      ray_dir_hor.obj   = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
      ray_point_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};

      // horizontal ray container
      ray_hor = {.direction = std::move(ray_dir_hor),    
                 .point     = std::move(ray_point_hor)};

      planes_hor.reserve(6* state.range(0)/Scalar_v::Size + 1);
      for (size_t s = 0; s < state.range(0)/Scalar_v::Size; s++) {
        pl_normals_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        pl_points_hor.obj  = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        planes_hor.push_back({.normals = std::move(pl_normals_hor), 
                              .points  = std::move(pl_points_hor)});
      }
      // padding at the end of data container needed (just add one more calculation for simplicity)
      if (state.range(0)/Scalar_v::Size != 0) {
        pl_normals_hor.obj = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        pl_points_hor.obj  = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
        planes_hor.push_back({.normals = std::move(pl_normals_hor), 
                              .points  = std::move(pl_points_hor)});
      }
    }

    void TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes_hor.clear();
      }
    }
};

//----------------------------------------------------Run Tests

//----------------------------------------------------//
// Eigen                                              //
//----------------------------------------------------//
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

//----------------------------------------------------//
// Vc                                                 //
//----------------------------------------------------//

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

//
// Use a horizontal vectorization and data set
//
BENCHMARK_DEFINE_F(HorizSetup, intersectVcHoriz)(benchmark::State& state){

  //TODO: Make pointer access threadsafe
  for (auto _: state) {
    #ifdef DEBUG 
    auto t1 = clock::now();
    #endif
    for (auto &plane : planes_hor) {
      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        vc_intersect_horiz<vector_v>(ray_hor, plane)
      );
      // Prevent compiler from optimizing away the loop
      benchmark::ClobberMemory();
    }

    #ifdef DEBUG 
    auto t2 = clock::now();
    auto duration = duration_cast<unit_ms>(t2 - t1);
    std::cout << "Eigen 4D (w vec): " << duration.count() << "ms\n";
    #endif
  }
  
}

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
    // Allow vector data to be clobbered in memory
    benchmark::DoNotOptimize(results.data());

    vc_intersect_horiz<vector_v>(ray_hor, planes_hor, results);

    // Prevent compiler from optimizing away the loop
    benchmark::ClobberMemory();
    //TODO: threadsafety
    results.clear();
  }
}

// Horizontal vectorized intersection is very sensitive to system load, so execute first
BENCHMARK_REGISTER_F(HorizSetup, intersectVcHoriz)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("VcHoriz")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef NO_MULTI_THREAD
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(HorizSetup, intersectVcHoriz_wres)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("VcHoriz_wres")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef NO_MULTI_THREAD
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(VertSetup, intersectEigen4D)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("Eigen4D")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef NO_MULTI_THREAD
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(VertSetup, intersectEigen4D_wres)
  ->DenseRange(surf_step, n_surf_steps*surf_step, surf_step)
  ->Name("Eigen4D_wres")
  //->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef NO_MULTI_THREAD
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

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


