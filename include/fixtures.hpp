#include "types.hpp"

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
constexpr size_t surf_step    = 10;
constexpr size_t n_surf_steps = 100;
constexpr size_t n_surf_mult = 5;
constexpr size_t n_surf_min = 10;
#ifdef NO_MULTI_THREAD
#define MULTI_THREAD 1
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

generator_t generator(42);

boost::uniform_real<Scalar> uni_dist(0,1000);
rand_t uni(generator, uni_dist);

//----------------------------------------------------Fill Data Types

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
    void SetUp(const ::benchmark::State& state);

    void TearDown(const ::benchmark::State& state);
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
    void SetUp(const ::benchmark::State& state);

    void TearDown(const ::benchmark::State& state);  
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
    void SetUp(const ::benchmark::State& state);

    void TearDown(const ::benchmark::State& state);
};

} // namespace g_bench

} //namespace vec_intr

#include <fixtures.ipp>