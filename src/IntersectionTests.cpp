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

//----------------------------------------------------Fill Data Types

  //using namespace vec_test;
struct data_setup {

  // Make sure the compiler does not optimize the computations away
  Scalar check_sum;

  Vector3<Scalar_v> ray_dir_hor, ray_point_hor;

  aligned::vector<vector_s>         pl_normals, pl_points;
  aligned::vector<Vector3<Scalar> > pl_points_struct, pl_normals_struct;
  aligned::vector<vector_v>         pl_normals_hor, pl_points_hor;

  ray_data<Vector4_s> ray;
  plane_data<aligned::vector<vector_s>> planes;

  ray_data<Vector3<Scalar_v>> ray_struct;
  plane_data<aligned::vector<Vector3<Scalar> >> planes_struct;

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

    // vertical data containers
    ray    = {.direction = ray_dir, 
              .point     = ray_point}; 
    planes = {.normals = pl_normals, 
              .points  = pl_points};

    // AoS container
    ray_struct    = {.direction = ray_dir_hor, 
                     .point     = ray_point_hor}; 
    planes_struct = {.normals = pl_normals_struct, 
                     .points  = pl_points_struct};

    // horizontal ray container
    ray_hor    = {.direction = ray_dir_hor,    
                  .point     = ray_point_hor};
    planes_hor = {.normals = pl_normals_hor, 
                  .points  = pl_points_hor};
  }
};
//----------------------------------------------------Run Tests

//----------------------------------------------------//
// Eigen                                              //
//----------------------------------------------------//
template <unsigned int kPlanes> void intersectEigen4D(ray_data<Vector4_s>& ray,
                                                      plane_data<aligned::vector<vector_s>>& planes) {
  Scalar check_sum = 0.0;
  auto padding = alignment % Scalar_v::Size;

  auto t1 = clock::now();
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
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Eigen 4D: " << duration.count() << "ms\n";
  std::cout << check_sum << std::endl;
}

template <unsigned int kPlanes> void intersectEigen4D_res(ray_data<Vector4_s>& ray,
                                                          plane_data<aligned::vector<vector_s>>& planes) {
  Scalar check_sum = 0.0;

  aligned::vector<intersection<Scalar, Vector4_s>> results;
  results.reserve(planes.points.size());

  auto t1 = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
      eig_intersect_4D<vector_s>(ray, planes, results);
      for (auto &intersection : results) check_sum += intersection.dist;
      #ifdef DEBUG
      if (nt % (tests-1)/2 == 0) {
        for (auto &intersection : results) {
          std::cout << "\n" << intersection.dist << std::endl;
          for (int i = 0; i < 4; i++) {
            std::cout << intersection.path[i] << ", ";
          }
        }
      }
      #endif
      results.clear();
  }
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Eigen 4D (w vec): " << duration.count() << "ms\n";
  std::cout << check_sum << std::endl;
}

template <unsigned int kPlanes> void intersectVcVert(ray_data<Vector4_s>& ray,
                                                      plane_data<aligned::vector<vector_s>>& planes) {
  Scalar check_sum = 0.0;
  
  auto t1 = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
    for (size_t i = 0; i < nSurfaces; i++) {
      auto intersection = vc_intersect_vert<Scalar_v, vector_s>(ray, planes.normals[i].obj, planes.points[i].obj);
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
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Vc vert: " << duration.count() << "ms\n";
  std::cout << check_sum << std::endl;
}

template <unsigned int kPlanes> void intersectVcVert_res(ray_data<Vector4_s>& ray,
                                                        plane_data<aligned::vector<vector_s>>& planes) {
  Scalar check_sum = 0.0;
  
  aligned::vector<intersection<Scalar, Vc::SimdArray<Scalar, 4>>> results;
  results.reserve(planes.points.size());

  auto t1 = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
    vc_intersect_vert<Scalar_v, vector_s>(ray, planes, results);
    for (auto &intersection : results) check_sum += intersection.dist;
    #ifdef DEBUG
    if (nt % (tests-1)/2 == 0) {
      for (auto &intersection : results) {
        std::cout << "\n" << intersection.dist << std::endl;
        for (int i = 0; i < 4; i++) {
          std::cout << intersection.path[i] << ", ";
        }
      }
    }
    #endif
    results.clear();
  }
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Vc vert (w vec): " << duration.count() << "ms\n";
  std::cout << check_sum << std::endl;
}


template <unsigned int kPlanes> void intersectVcHybrid(ray_data<Vector3<Scalar_v>>& ray,
                                                       plane_data<aligned::vector<Vector3<Scalar>>>& planes) {
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  auto t1 = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
    for (Index_v i(Scalar_v::IndexesFromZero()); (i < Index_v(planes.points.size())).isFull(); i += Index_v(Scalar_v::Size)) {
      Scalar_v pns_x = planes.normals[i][&Vector3<Scalar>::x];
      Scalar_v pns_y = planes.normals[i][&Vector3<Scalar>::y];
      Scalar_v pns_z = planes.normals[i][&Vector3<Scalar>::z];
      Vector3<Scalar_v> pl_normal_strc {.x = pns_x, .y = pns_y, .z = pns_z};

      Scalar_v pps_x = planes.points[i][&Vector3<Scalar>::x];
      Scalar_v pps_y = planes.points[i][&Vector3<Scalar>::y];
      Scalar_v pps_z = planes.points[i][&Vector3<Scalar>::z];
      Vector3<Scalar_v> pl_point_strc {.x = pps_x, .y = pps_y, .z = pps_z};
      plane_data<Vector3<Scalar_v>> planes = {.normals = pl_normal_strc, .points = pl_point_strc};

      auto intersection = vc_intersect_hybrid<Scalar_v>(ray, planes);
      check_sum_v += intersection.dist;
      #ifdef DEBUG
      if (nt % (tests-1)/2 == 0) {
        std::cout << intersection.dist << std::endl;
        std::cout << intersection.path.x << "\t" << intersection.path.y << "\t" << intersection.path.z << std::endl;
      }
      #endif
    }
  }
  check_sum = check_sum_v.sum();
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Vc hybrd: " << duration.count() << "ms\n";
  std::cout << check_sum << std::endl;
}

template <unsigned int kPlanes> void intersectVcHybrid_res(ray_data<Vector3<Scalar_v>>& ray,
                                                           plane_data<aligned::vector<Vector3<Scalar>>>& planes) {
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  aligned::vector<intersection<Scalar_v, Vector3<Scalar_v>>> results;
  results.reserve(planes.normals.size());
  
  auto t1 = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
      vc_intersect_hybrid<Scalar_v>(ray, planes, results);
      for (auto &intersection : results) check_sum_v += intersection.dist;
      #ifdef DEBUG
      if (nt % (tests-1)/2 == 0) {
        for (auto &intersection : results) {
          std::cout << intersection.dist << std::endl;
          std::cout << intersection.path.x << "\t" << intersection.path.y << "\t" << intersection.path.z << std::endl;
        }
      }
      #endif
      results.clear();
  }
  check_sum = check_sum_v.sum();
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Vc hybrd (w vec): " << duration.count() << "ms\n";
  std::cout << check_sum << std::endl;
}


template <unsigned int kPlanes> void intersectVcHoriz(ray_data<Vector3<Scalar_v>>& ray,
                                                      plane_data<aligned::vector<vector_v>>& planes) {
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  auto padding = planes.points.front().padding();

  // Access to raw data that will be loaded as scalar_v
  size_t n_bytes = planes.points.size() * (planes.points.front().n_elemts() + padding);
  size_t offset  = Scalar_v::Size + padding;
  if (n_bytes % (3*offset) != 0) std::cout << "Warning: Input container size is not a multiple simd vector size." << std::endl;
  size_t n_vec = n_bytes / (3*offset);

  auto t1 = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
      auto pl_normals_ptr = const_cast<const Scalar*>(planes.normals.front().data());
      auto pl_points_ptr  = const_cast<const Scalar*>(planes.points.front().data());

      for (size_t i = 0; i < n_vec; i++) {
        auto intersection = vc_intersect_horiz<Scalar_v, const Scalar*>(ray, pl_normals_ptr, pl_points_ptr, offset);

        check_sum_v     += intersection.dist;
        pl_normals_ptr += 3 * offset;
        pl_points_ptr  += 3 * offset;

        #ifdef DEBUG 
        if (nt % (tests-1)/2 == 0) {
          std::cout << intersection.dist << std::endl;
          std::cout << intersection.path.x << "\t" << intersection.path.y << "\t"
                    << intersection.path.z << std::endl;
        }
        #endif
      }
  }
  check_sum = check_sum_v.sum();
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Vc horizontal: " << duration.count() << "ms" << std::endl;
  std::cout << check_sum << std::endl;
}

template <unsigned int kPlanes> void intersectVcHoriz_res(ray_data<Vector3<Scalar_v>>& ray,
                                                          plane_data<aligned::vector<vector_v>>& planes) {
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  auto padding = planes.points.front().padding();

  aligned::vector<intersection<Scalar_v, Vector3<Scalar_v>>> results;
  results.reserve(planes.points.size() * planes.points.front().n_elemts()/Scalar_v::Size);

  auto t1 = clock::now();
  for (size_t nt = 0; nt < tests; ++nt) {
      vc_intersect_horiz<Scalar_v, vector_v::obj_type>(ray, planes, results, padding);
      for (auto &intersection : results) check_sum_v += intersection.dist;
      #ifdef DEBUG 
      if (nt % (tests-1)/2 == 0) {
        for (auto &intersection : results) {
          std::cout << intersection.dist << std::endl;
          std::cout << intersection.path.x << "\t" << intersection.path.y << "\t"
                    << intersection.path.z << std::endl;
        }
      }
      #endif
      results.clear();
  }
  check_sum = check_sum_v.sum();
  auto t2 = clock::now();
  auto duration = duration_cast<unit_ms>(t2 - t1);
  std::cout << "Vc horizontal (w vec): " << duration.count() << "ms" << std::endl;
  std::cout << check_sum << std::endl;
}


//----------------------------------------------------Run Boost Tests

BOOST_FIXTURE_TEST_SUITE(VectIntersect, data_setup)


BOOST_AUTO_TEST_CASE(IntersectEigen4D)    {intersectEigen4D<nSurfaces>(ray, planes);}
BOOST_AUTO_TEST_CASE(IntersectEigen4D_res){intersectEigen4D_res<nSurfaces>(ray, planes);}


BOOST_AUTO_TEST_CASE(IntersectVcVert)     {intersectVcVert<nSurfaces>(ray, planes);}
BOOST_AUTO_TEST_CASE(IntersectVcVert_res) {intersectVcVert_res<nSurfaces>(ray, planes);}


BOOST_AUTO_TEST_CASE(IntersectVcHybrid)     {intersectVcHybrid<nSurfaces>(ray_struct, 
                                                                          planes_struct);}
BOOST_AUTO_TEST_CASE(IntersectVcHybrid_res) {intersectVcHybrid_res<nSurfaces>(ray_struct, 
                                                                              planes_struct);}


BOOST_AUTO_TEST_CASE(IntersectVcHoriz)    {intersectVcHoriz<nSurfaces>(ray_hor, 
                                                                       planes_hor);}
BOOST_AUTO_TEST_CASE(IntersectVcHorizres) {intersectVcHoriz_res<nSurfaces>(ray_hor, 
                                                                           planes_hor);}


BOOST_AUTO_TEST_SUITE_END()

} // namespace test

} //namespace vec_intr


