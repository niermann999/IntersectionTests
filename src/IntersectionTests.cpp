 /**
 * author: asalzburger@gmail.com
 **/
#include <types.hpp>
#include <intersectors.hpp>

#include <array>
#include <ctime>            // std::time
#include <chrono>
#include <iostream>
#include <limits>

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

namespace vec_intr {

namespace test {


// This is a typedef for a random number generator.
// Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
using generator_t = boost::minstd_rand;
// Define a uniform random number distribution which produces "double"
// values between 0 and 1 (0 inclusive, 1 exclusive).
using rand_t = boost::variate_generator<generator_t&, boost::uniform_real<Scalar> >;

#ifdef DEBUG
constexpr size_t tests = 100;
constexpr size_t nSurfaces = 8;
#else
constexpr size_t tests = 100;
constexpr size_t nSurfaces = 128000;
#endif

// Make sure the memory layout is compatible with Vc Vectors and set corresponding LA wrappers as types
static_assert(data_trait<Vector4_s>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");
static_assert(data_trait<VectorV_s>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");
//static_assert(data_trait<Transform4>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");

using vector_s = data_trait<Vector4_s>::type;
using vector_v = data_trait<VectorV_s>::type;
//using transf_s = data_trait<Transform4>::type;

// Measure execution time
using clock = std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using unit_ms = std::chrono::milliseconds;

 // Allow Vc to get alignment right
  //template <typename T, typename Allocator = std::allocator<T>>
  //template <typename T, typename Allocator = Eigen::aligned_allocator<T>>
  //template <typename T, typename Allocator = Vc::Allocator<T>>
  // Add subscript operator to allow for gather operations from AoS
  //using vector_a = Vc::Common::AdaptSubscriptOperator<std::vector<T, Allocator>>;

//----------------------------------------------------Fill Data Types

  //using namespace vec_test;
struct data_setup {

  // Make sure the compiler does not optimize the computations away
  Scalar check_sum;

  Vector3<Scalar_v> ray_v, ray_p;
  Vector3<Scalar_v> ray_dir_hor, ray_point_hor;

  aligned::vector<Vector3<Scalar> > pl_points_struct, pl_normals_struct;
  aligned::vector<vector_v>         pl_normals_hor, pl_points_hor;

  ray_data<Vector4_s> ray;
  plane_data<aligned::vector<vector_s>> planes;

  ray_data<Vector3<Scalar_v>> ray_hor;
  plane_data<aligned::vector<vector_v>> planes_hor;

  data_setup () {
    // Produce test data (either predefined or random)
    #ifdef DEBUG

    // Ray to be intersected with different plains
    Scalar data[4] = {0.0, -1.0, -1.0, 0.0};//, 0.0, 0.0, 0.0, 0.0, 0.0};
    Vector4_s ray_dir = Vector4_s(data);
    data[2] = 10.0;
    data[1] = 0.0;
    Vector4_s ray_point = Vector4_s(data);
    data[2] = 1.0;

    // Depends on scalar precision and vec standard (avx...)
    /*Scalar r_data[8] = {0.0, -1.0, -1.0, 0.0, 0.0, -1.0, -1.0, 0.0};
    // Same starting position
    VectorV_s ray_dir = VectorV_s(r_data);
    r_data[2] = 10.0;
    r_data[1] = 0.0;
    r_data[6] = 10.0;
    r_data[5] = 0.0;
    data[2] = 10.0;
    data[1] = 0.0;
    VectorV_s rp = VectorV_s(r_data);*/
    // Same starting position

    // For the moment same normal vectors
    vector_s pl_normal {.obj = vector_s::obj_type(data)};

    // plane normals and points
    data[2] = 5.0;
    aligned::vector<vector_s> pl_points;
    aligned::vector<vector_s> pl_normals;
    pl_points.reserve(nSurfaces);
    pl_normals.reserve(nSurfaces);
    for (size_t i = 0; i < nSurfaces; i++) {
      pl_points.push_back({.obj = vector_s::obj_type(data)});
      pl_normals.push_back(pl_normal);
      data[2]++;
    }

    // AoS data
    ray_v = {.x=Scalar_v(0.0), .y=Scalar_v(-1.0), .z=Scalar_v(-1.0)};
    ray_p = {.x=Scalar_v(0.0), .y=Scalar_v(0.0),  .z=Scalar_v(10.0)};

    Vector3<Scalar> pp0 {.x=0.0, .y=0.0, .z=5.0};
    Vector3<Scalar> pp1 {.x=0.0, .y=0.0, .z=6.0};
    Vector3<Scalar> pp2 {.x=0.0, .y=0.0, .z=7.0};
    Vector3<Scalar> pp3 {.x=0.0, .y=0.0, .z=8.0};
    Vector3<Scalar> pp4 {.x=0.0, .y=0.0, .z=9.0};
    Vector3<Scalar> pp5 {.x=0.0, .y=0.0, .z=10.0};
    Vector3<Scalar> pp6 {.x=0.0, .y=0.0, .z=11.0};
    Vector3<Scalar> pp7 {.x=0.0, .y=0.0, .z=12.0};
    Vector3<Scalar> pp8 {.x=0.0, .y=0.0, .z=13.0};

    Vector3<Scalar> plain_n {.x=0.0, .y=0.0, .z=1.0};

    pl_points_struct = {pp0, pp1, pp2, pp3, pp4, pp5, pp6, pp7};
    pl_normals_struct = {plain_n, plain_n, plain_n, plain_n, 
                                                    plain_n, plain_n, plain_n, plain_n};

    // horizontal data (interleaved)
    ray_dir_hor   = {.x= Scalar_v(0.0), .y=Scalar_v(-1.0), .z=Scalar_v(-1.0)};
    ray_point_hor = {.x= Scalar_v(0.0), .y=Scalar_v(0.0), .z=Scalar_v(10.0)};


    for (size_t offset = 0; offset < nSurfaces/Scalar_v::Size; offset++) {
      for (size_t i = 0; i < 3; i++) {
        vector_v pl_points_i {};
        vector_v pl_normals_i {};
        for (size_t j = 0; j < Scalar_v::Size; j++) {
          pl_points_i.obj[j] = pl_points[j+offset*Scalar_v::Size].obj[i];
          pl_normals_i.obj[j] = pl_normals[j+offset*Scalar_v::Size].obj[i];
        }
        pl_points_hor.push_back(pl_points_i);
        pl_normals_hor.push_back(pl_normals_i);
      }
    }
    

    for (const auto& vec : pl_points_hor) {
      std::cout << vec.obj.transpose() <<" \t";
    }
    std::cout << std::endl;
    for (const auto& vec : pl_normals_hor) {
      std::cout << vec.obj.transpose() << " \t";
    }
    std::cout << std::endl;
    
    #else

    //-----------------------------------------------Random data
    
    generator_t generator(42);

    boost::uniform_real<Scalar> uni_dist(0,1000);
    rand_t uni(generator, uni_dist);

    // Same starting position
    Vector4_s ray_dir = Vector4_s::Random();
    Vector4_s ray_point = Vector4_s::Random();

    aligned::vector<vector_s> pl_normals;
    aligned::vector<vector_s> pl_points;
    pl_normals.reserve(nSurfaces);
    pl_points.reserve(nSurfaces);
    for (size_t i = 0; i < nSurfaces; i++) {
      pl_normals.push_back({.obj = vector_s::obj_type::Random()});
      pl_points.push_back({.obj = vector_s::obj_type::Random()});
    }

    // AoS data
    aligned::vector<Vector3<Scalar> > pl_normals_struct;
    aligned::vector<Vector3<Scalar> > pl_points_struct;
    pl_normals_struct.reserve(nSurfaces);
    pl_points_struct.reserve(nSurfaces);

    for (size_t i = 0; i < nSurfaces; i++) {
      pl_normals_struct.push_back({.x=uni(), .y=uni(), .z=uni()});
      pl_points_struct.push_back({.x=uni(), .y=uni(), .z=uni()});
    }

    // Horizontal data (interleaved)
    Vector3<Scalar_v> ray_dir_hor   {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};
    Vector3<Scalar_v> ray_point_hor {.x= Scalar_v(uni()), .y=Scalar_v(uni()), .z=Scalar_v(uni())};

    aligned::vector<vector_v> pl_normals_hor;
    aligned::vector<vector_v> pl_points_hor;
    // dimension * number of matrices needed
    for (size_t offset = 0; offset < 3 * nSurfaces/Scalar_v::Size; offset++) {
      pl_normals_hor.push_back({.obj = vector_v::obj_type::Random()});
      pl_points_hor.push_back({.obj = vector_v::obj_type::Random()});
    }
    
    /*for (const auto& vec : pl_points_hor) {
      std::cout << vec.obj.transpose() <<" \t";
    }
    std::cout << std::endl;
    for (const auto& vec : pl_normals_hor) {
      std::cout << vec.obj.transpose() << " \t";
    }
    std::cout << std::endl;*/

    #endif

    // initialize vertical data containers
    ray = {.direction = ray_dir, .point     = ray_point}; 

    planes = {.normals = pl_normals, .points  = pl_points};

    // horizontal ray container
    ray_data<Vector3<Scalar_v>> ray_hor = {.direction = ray_dir_hor,
                                          .point     = ray_point_hor};

    plane_data<aligned::vector<vector_v>> planes_hor = {.normals = pl_normals_hor,
                                                        .points  = pl_points_hor};
  }
};
//----------------------------------------------------Run Tests

template <unsigned int kPlanes> void intersectEigen4D(ray_data<Vector4_s> ray,
                                                      plane_data<aligned::vector<vector_s>> planes,
                                                      Scalar check_sum) {
  // Just the intersections
  auto padding = alignment % Scalar_v::Size;

  auto t1_eig = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
    for (size_t i = 0; i < nSurfaces; i++) {
      auto intersection = eig_intersect_4D(ray, planes.normals[i].obj, planes.points[i].obj);
      check_sum += intersection.dist;
      #ifdef DEBUG
      if (nt % (tests-1)/2 == 0) {
        std::cout << "\n" << intersection.dist << std::endl;
        for (int i = 0; i < 4; i++) {
          std::cout << intersection.path[i] << ", ";
        }
      }
      #endif
    }
  }
  auto t2_eig = clock::now();
  auto duration_eig = duration_cast<unit_ms>(t2_eig - t1_eig);
  std::cout << "Eigen 4D: " << duration_eig.count() << "ms\n";
  std::cout << check_sum << std::endl;
}

template <unsigned int kPlanes> void intersectEigen4D_res(ray_data<Vector4_s> ray,
                                                          plane_data<aligned::vector<vector_s>> planes,
                                                          Scalar check_sum) {
  aligned::vector<intersection<Scalar, Vector4_s>> results_scalar;
  results_scalar.reserve(planes.points.size());

  auto t1_eig_wres = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
      eig_intersect_4D<vector_s>(ray, planes, results_scalar);
      for (auto &intersection : results_scalar) check_sum += intersection.dist;
      #ifdef DEBUG
      if (nt % (tests-1)/2 == 0) {
        for (auto &intersection : results_scalar) {
          std::cout << "\n" << intersection.dist << std::endl;
          for (int i = 0; i < 4; i++) {
            std::cout << intersection.path[i] << ", ";
          }
        }
      }
      #endif
      results_scalar.clear();
  }
  auto t2_eig_wres = clock::now();
  auto duration_eig_wres = duration_cast<unit_ms>(t2_eig_wres - t1_eig_wres);
  std::cout << "Eigen 4D (w vec): " << duration_eig_wres.count() << "ms\n";
  std::cout << check_sum << std::endl;
}


//----------------------------------------------------Run Boost Tests

BOOST_FIXTURE_TEST_SUITE(VectIntersect, data_setup)


BOOST_AUTO_TEST_CASE(SingleIntersection4D) { intersectEigen4D<nSurfaces>(ray, planes, check_sum); }
BOOST_AUTO_TEST_CASE(SingleIntersection4D_res) { intersectEigen4D_res<nSurfaces>(ray, planes, check_sum); }


BOOST_AUTO_TEST_SUITE_END()

} // namespace test

} //namespace vec_intr


