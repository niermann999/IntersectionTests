 /**
 * author: joana.niermann@cern.ch
 **/
#include "types.hpp"
#include "intersectors.hpp"

#ifdef DEBUG
#include <chrono>
#include <ctime>
#endif
#include <atomic>
//#include <barrier> //c++ 20
#include <iostream>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

/*#include <boost/thread.hpp>
#include <boost/thread/barrier.hpp>
#include <boost/bind.hpp>
#include <boost/atomic.hpp>*/

#include <benchmark/benchmark.h>

namespace vec_intr {

namespace g_bench {

// Iterations of a benchmark
constexpr size_t gbench_test_itrs = 10000;
// Repetitions of a benchmark
constexpr size_t gbench_test_repts = 10;
// Number of rand. gen. surfaces to intersect
constexpr size_t nSurfaces = 132;
#ifdef DEBUG
constexpr size_t nThreads  = 1;
#endif

// Make sure the memory layout is compatible with Vc Vectors and set corresponding LA wrappers as types
static_assert(data_trait<Vector4_s>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");
static_assert(data_trait<VectorV_s>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");
//static_assert(data_trait<Transform4>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");

using vector_s = data_trait<Vector4_s>::type;
using vector_v = data_trait<VectorV_s>::type;
//using transf_s = data_trait<Transform4>::type;


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
    Scalar   check_sum;

    aligned::vector<vector_s> pl_normals, pl_points;
    ray_data<Vector4_s> ray;
    plane_data<aligned::vector<vector_s>> planes;

    Vector4_s ray_dir, ray_point;

    #ifdef DEBUG
    // manual time measurement using chrono
    unit_ms duration;
    #endif

    ~VertSetup() = default;
    void SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      ray_dir   = Vector4_s::Random();
      ray_point = Vector4_s::Random();

      pl_normals.reserve(state.range(0));
      pl_points.reserve(state.range(0));
      for (size_t i = 0; i < state.range(0); i++) {
        pl_normals.push_back({.obj = vector_s::obj_type::Random()});
        pl_points.push_back( {.obj = vector_s::obj_type::Random()});
      }

        // vertical data containers
      ray    = {.direction = std::move(ray_dir), 
                .point     = std::move(ray_point)}; 
      planes = {.normals = std::move(pl_normals), 
                .points  = std::move(pl_points)};
    }

    void TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes.normals.clear();
        planes.points.clear();
      }
    }
    
};

// Vertical data contained in structs
class HybridSetup : public benchmark::Fixture {

    public:
    Scalar   check_sum;
    Scalar_v check_sum_v;

    Vector3<Scalar_v> ray_dir_hor, ray_point_hor;
    aligned::vector<Vector3<Scalar> > pl_points_struct, pl_normals_struct;

    ray_data<Vector3<Scalar_v>> ray_struct;
    plane_data<aligned::vector<Vector3<Scalar> >> planes_struct;

    ray_data<Vector3<Scalar_v>> ray_hor;
    plane_data<aligned::vector<vector_v>> planes_hor;

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
      ray_dir_hor   = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
      ray_point_hor = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};

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
    Scalar   check_sum;
    Scalar_v check_sum_v;

    Vector3<Scalar_v>         ray_dir_hor, ray_point_hor;
    aligned::vector<vector_v> pl_normals_hor, pl_points_hor;

    ray_data<Vector3<Scalar_v>> ray_hor;
    plane_data<aligned::vector<vector_v>> planes_hor;

    #ifdef DEBUG
    // manual time measurement using chrono
    unit_ms duration;
    #endif

    ~HorizSetup() = default;
    void SetUp(const ::benchmark::State& state) {
      // Only the first thread sets the test env up
      if (state.thread_index != 0) return;

      // Horizontal data (interleaved)
      ray_dir_hor  = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
      ray_point_hor = {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};

      // need 3 vc vectors every "vector-reg. width" of surfaces (can compute "vector-reg. width" surfaces at the same time)
      pl_normals_hor.reserve(3* state.range(0)/Scalar_v::Size + 1);
      pl_points_hor.reserve(3* state.range(0)/Scalar_v::Size + 1);
      for (size_t s = 0; s < state.range(0)/Scalar_v::Size; s++) {
        pl_normals_hor.push_back({.obj = vector_v::obj_type::Random()}); //x
        pl_normals_hor.push_back({.obj = vector_v::obj_type::Random()}); //y
        pl_normals_hor.push_back({.obj = vector_v::obj_type::Random()}); //z

        pl_points_hor.push_back({.obj = vector_v::obj_type::Random()});
        pl_points_hor.push_back({.obj = vector_v::obj_type::Random()});
        pl_points_hor.push_back({.obj = vector_v::obj_type::Random()});
      }
      // padding at the end of data container needed (just add one more calculation for simplicity)
      if (state.range(0)/Scalar_v::Size != 0) {
        pl_normals_hor.push_back({.obj = vector_v::obj_type::Random()}); //x
        pl_normals_hor.push_back({.obj = vector_v::obj_type::Random()}); //y
        pl_normals_hor.push_back({.obj = vector_v::obj_type::Random()}); //z

        pl_points_hor.push_back({.obj = vector_v::obj_type::Random()});
        pl_points_hor.push_back({.obj = vector_v::obj_type::Random()});
        pl_points_hor.push_back({.obj = vector_v::obj_type::Random()});
      }

      // horizontal ray container
      ray_hor    = {.direction = std::move(ray_dir_hor),    
                    .point     = std::move(ray_point_hor)};
      planes_hor = {.normals = std::move(pl_normals_hor), 
                    .points  = std::move(pl_points_hor)};
    }

    void TearDown(const ::benchmark::State& state) {
      // Only one thread frees recources
      if (state.thread_index == 0) {
        planes_hor.normals.clear();
        planes_hor.points.clear();
      }
    }
};

//----------------------------------------------------Run Tests

//----------------------------------------------------//
// Eigen                                              //
//----------------------------------------------------//
BENCHMARK_DEFINE_F(VertSetup, intersectEigen4D)(benchmark::State& state) {

  for (auto _: state) {
    for (size_t i = 0; i < state.range(0); i++) {
      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        eig_intersect_4D(ray, planes.normals[i].obj, planes.points[i].obj)
      );
      // Prevent compiler from optimizing away the loop
      benchmark::ClobberMemory();
    }
  }
}

BENCHMARK_DEFINE_F(VertSetup, intersectEigen4D_wres)(benchmark::State& state) {
  aligned::vector<intersection<Scalar, Vector4_s>> results;
  if (state.thread_index == 0) results.reserve(planes.points.size());

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
    for (size_t i = 0; i < state.range(0); i++) {
      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        vc_intersect_vert<Scalar_v, vector_s>(ray, planes.normals[i].obj, planes.points[i].obj)
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

  aligned::vector<intersection<Scalar, Vc::SimdArray<Scalar, 4>>> results;
  if (state.thread_index == 0) results.reserve(planes.points.size());

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    // Allow vector data to be clobbered in memory
    benchmark::DoNotOptimize(results.data());

    vc_intersect_vert<Scalar_v, vector_s>(ray, planes, results);

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

  for (auto _: state) {
    for (Index_v i(Scalar_v::IndexesFromZero()); (i < Index_v(planes_struct.points.size())).isFull(); i += Index_v(Scalar_v::Size)) {
      Scalar_v pns_x = planes_struct.normals[i][&Vector3<Scalar>::x];
      Scalar_v pns_y = planes_struct.normals[i][&Vector3<Scalar>::y];
      Scalar_v pns_z = planes_struct.normals[i][&Vector3<Scalar>::z];
      Vector3<Scalar_v> pl_normal_strc {.x = pns_x, .y = pns_y, .z = pns_z};

      Scalar_v pps_x = planes_struct.points[i][&Vector3<Scalar>::x];
      Scalar_v pps_y = planes_struct.points[i][&Vector3<Scalar>::y];
      Scalar_v pps_z = planes_struct.points[i][&Vector3<Scalar>::z];
      Vector3<Scalar_v> pl_point_strc {.x = pps_x, .y = pps_y, .z = pps_z};
      plane_data<Vector3<Scalar_v>> planes_strcts = {.normals = pl_normal_strc, .points = pl_point_strc};

      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        vc_intersect_hybrid<Scalar_v>(ray_struct, planes_strcts)
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
  
  aligned::vector<intersection<Scalar_v, Vector3<Scalar_v>>> results;
  if (state.thread_index == 0)  results.reserve(planes_struct.normals.size());

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    // Allow vector data to be clobbered in memory
    benchmark::DoNotOptimize(results.data());

    vc_intersect_hybrid<Scalar_v>(ray_struct, planes_struct, results);

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
    size_t padding = planes_hor.points.front().padding();
    // Access to raw data that will be loaded as scalar_v
    size_t n_float_pnt = planes_hor.points.size() * (planes_hor.points.front().n_elemts() + padding);
    size_t offset  = planes_hor.points.front().n_elemts() + padding;
    size_t n_loops = n_float_pnt / (3*offset);
    #ifdef DEBUG
    // Process 3 geometrical coordinates
    if (n_float_pnt % (3*offset) != 0) std::cout << "Warning: Input container size is not a multiple simd vector size." << std::endl;
    #endif

    auto pl_normals_ptr = const_cast<const Scalar*>(planes_hor.normals.front().data());
    auto pl_points_ptr  = const_cast<const Scalar*>(planes_hor.points.front().data());

    #ifdef DEBUG 
    auto t1 = clock::now();
    #endif
    for (size_t i = 0; i < n_loops; i++) {
      // Allow return value to be clobbered in memory
      benchmark::DoNotOptimize(
        vc_intersect_horiz<Scalar_v, const Scalar*>(ray_hor, pl_normals_ptr, pl_points_ptr, offset)
      );
      // Prevent compiler from optimizing away the loop
      benchmark::ClobberMemory();

      pl_normals_ptr += 3 * offset;
      pl_normals_ptr += 3 * offset;
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
  aligned::vector<intersection<Scalar_v, Vector3<Scalar_v>>> results;
  if (state.thread_index == 0) {
    results.reserve(planes_hor.points.size() * planes_hor.points.front().n_elemts()/Scalar_v::Size);
  }

  //TODO: Make vector operations threadsafe
  for (auto _: state) {
    size_t padding = planes_hor.points.front().padding();
    // Allow vector data to be clobbered in memory
    benchmark::DoNotOptimize(results.data());

    vc_intersect_horiz<Scalar_v, vector_v::obj_type>(ray_hor, planes_hor, results, padding);

    // Prevent compiler from optimizing away the loop
    benchmark::ClobberMemory();
    //TODO: threadsafety
    results.clear();
  }
}

BENCHMARK_REGISTER_F(VertSetup, intersectEigen4D)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("Eigen4D")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(VertSetup, intersectEigen4D_wres)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("Eigen4D_wres")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(VertSetup, intersectVcVert)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("VcVert")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(VertSetup, intersectVcVert_wres)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("VcVert_wres")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(HybridSetup, intersectVcHybrid)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("VcHybrid")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(HybridSetup, intersectVcHybrid_wres)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("VcHybrid_wres")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(HorizSetup, intersectVcHoriz)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("VcHoriz")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_REGISTER_F(HorizSetup, intersectVcHoriz_wres)
  ->DenseRange(0, 10*nSurfaces, nSurfaces)
  ->Name("VcHoriz_wres")
  ->Iterations(gbench_test_itrs)
  ->Repetitions(gbench_test_repts)
  ->DisplayAggregatesOnly(true)
  #ifdef DEBUG
  ->Threads(nThreads);
  #else
  ->ThreadPerCpu();
  #endif

BENCHMARK_MAIN();

} // namespace g_bench

} //namespace vec_intr


