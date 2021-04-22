 /**
 * author: joana.niermann@cern.ch
 **/
#include <types.hpp>
#include <intersectors.hpp>

#include <iostream>
#include <limits>

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

namespace vec_intr {

namespace test {

constexpr size_t nTests = 100;
constexpr size_t nSurfaces = 8;

// Make sure the memory layout is compatible with Vc Vectors and set corresponding LA wrappers as types
static_assert(data_trait<Vector4_s>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");
static_assert(data_trait<VectorV_s>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");
//static_assert(data_trait<Transform4>::is_vec_layout, "Input type has non-compatible memory layout for vectorization");

using vector_s = data_trait<Vector4_s>::type;
using vector_v = data_trait<VectorV_s>::type;
//using transf_s = data_trait<Transform4>::type;

//----------------------------------------------------Fill Data Types

// Produce test data
struct data_setup {
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
    // Ray to be intersected with different plains
    Scalar data[4] = {0.0, -1.0, -1.0, 0.0};//, 0.0, 0.0, 0.0, 0.0, 0.0};
    Vector4_s ray_dir = Vector4_s(data);
    data[2] = 10.0;
    data[1] = 0.0;
    Vector4_s ray_point = Vector4_s(data);
    data[2] = 1.0;

    // For the moment same normal vectors
    vector_s pl_normal {.obj = vector_s::obj_type(data)};

    // plane normals and points
    data[2] = 5.0;
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
  std::cout << "Running Eigen 4D test ...";
  Scalar check_sum = 0.0;

  for (size_t nt = 0; nt < nTests; ++nt) {
    for (size_t i = 0; i < nSurfaces; i++) {
      auto intersection = eig_intersect_4D(ray, planes.normals[i].obj, planes.points[i].obj);
      check_sum += intersection.dist;
      #ifdef DEBUG
      if (nt % (nTests-1)/2 == 0) {
        std::cout << "\n" << intersection.dist << std::endl;
        for (int i = 0; i < 4; i++) {
          std::cout << intersection.path[i] << ", ";
        }
      }
      #endif
    }
  }
  // TODO watch the floating point comparision!!
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\t\t\t\t\tdone" << std::endl;
}

template <unsigned int kPlanes> void intersectEigen4D_res(ray_data<Vector4_s>& ray,
                                                          plane_data<aligned::vector<vector_s>>& planes) {
  std::cout << "Running Eigen with vector container 4D test ...";
  Scalar check_sum = 0.0;

  aligned::vector<intersection<Scalar, Vector4_s>> results;
  results.reserve(planes.points.size());

  for (size_t nt = 0; nt < nTests; ++nt) {
      eig_intersect_4D<vector_s>(ray, planes, results);
      for (auto &intersection : results) check_sum += intersection.dist;
      #ifdef DEBUG
      if (nt % (nTests-1)/2 == 0) {
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
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\t\t\tdone" << std::endl;
}

//----------------------------------------------------//
// Vc                                                 //
//----------------------------------------------------//

//
// Use simdArray on unchanged data set
//
template <unsigned int kPlanes> void intersectVcVert(ray_data<Vector4_s>& ray,
                                                      plane_data<aligned::vector<vector_s>>& planes) {
  std::cout << "Running vert. Vectoriztion test ...";
  Scalar check_sum = 0.0;
  
  for (size_t nt = 0; nt < nTests; ++nt) {
    for (size_t i = 0; i < nSurfaces; i++) {
      auto intersection = vc_intersect_vert<Scalar_v, vector_s>(ray, planes.normals[i].obj, planes.points[i].obj);
      check_sum += intersection.dist;
      #ifdef DEBUG
      if (nt % (nTests-1)/2 == 0) {
        std::cout << "\n" << intersection.dist << std::endl;
        for (int i = 0; i < 4; i++) {
          std::cout << intersection.path[i] << ", ";
        }
      }
      #endif
    }
  }
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\t\t\t\tdone" << std::endl;
}

//
// Use simdArray on unchanged data set and save the results in a container
//
template <unsigned int kPlanes> void intersectVcVert_res(ray_data<Vector4_s>& ray,
                                                        plane_data<aligned::vector<vector_s>>& planes) {
  std::cout << "Running vert. Vectoriztion with vector container test ...";
  Scalar check_sum = 0.0;
  
  aligned::vector<intersection<Scalar, Vc::SimdArray<Scalar, 4>>> results;
  results.reserve(planes.points.size());

  for (size_t nt = 0; nt < nTests; ++nt) {
    vc_intersect_vert<Scalar_v, vector_s>(ray, planes, results);
    for (auto &intersection : results) check_sum += intersection.dist;
    #ifdef DEBUG
    if (nt % (nTests-1)/2 == 0) {
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
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\tdone" << std::endl;
}

//
// Use a gather on a structured data set then vectorize horizontaly
//
template <unsigned int kPlanes> void intersectVcHybrid(ray_data<Vector3<Scalar_v>>& ray,
                                                       plane_data<aligned::vector<Vector3<Scalar>>>& planes) {
  std::cout << "Running hybrid Vectoriztion test ...";
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  for (size_t nt = 0; nt < nTests; ++nt) {
    for (Index_v i(Scalar_v::IndexesFromZero()); (i < Index_v(planes.points.size())).isFull(); i += Index_v(Scalar_v::Size)) {
      Scalar_v pns_x = planes.normals[i][&Vector3<Scalar>::x];
      Scalar_v pns_y = planes.normals[i][&Vector3<Scalar>::y];
      Scalar_v pns_z = planes.normals[i][&Vector3<Scalar>::z];
      Vector3<Scalar_v> pl_normal_strc {.x = pns_x, .y = pns_y, .z = pns_z};

      Scalar_v pps_x = planes.points[i][&Vector3<Scalar>::x];
      Scalar_v pps_y = planes.points[i][&Vector3<Scalar>::y];
      Scalar_v pps_z = planes.points[i][&Vector3<Scalar>::z];
      Vector3<Scalar_v> pl_point_strc {.x = pps_x, .y = pps_y, .z = pps_z};
      plane_data<Vector3<Scalar_v>> planes_strcts = {.normals = pl_normal_strc, .points = pl_point_strc};

      auto intersection = vc_intersect_hybrid<Scalar_v>(ray, planes_strcts);
      check_sum_v += intersection.dist;
      #ifdef DEBUG
      if (nt % (nTests-1)/2 == 0) {
        std::cout << intersection.dist << std::endl;
        std::cout << intersection.path.x << "\t" << intersection.path.y << "\t" << intersection.path.z << std::endl;
      }
      #endif
    }
  }
  check_sum = check_sum_v.sum();
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\t\t\t\tdone" << std::endl;
}

//
// Use a gather on a structured data set then vectorize horizontaly 
// and save the results in a container
//
template <unsigned int kPlanes> void intersectVcHybrid_res(ray_data<Vector3<Scalar_v>>& ray,
                                                           plane_data<aligned::vector<Vector3<Scalar>>>& planes) {
  std::cout << "Running hybrid Vectoriztion with vector container test ...";
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  aligned::vector<intersection<Scalar_v, Vector3<Scalar_v>>> results;
  results.reserve(planes.normals.size());
  
  for (size_t nt = 0; nt < nTests; ++nt) {
      vc_intersect_hybrid<Scalar_v>(ray, planes, results);
      for (auto &intersection : results) check_sum_v += intersection.dist;
      #ifdef DEBUG
      if (nt % (nTests-1)/2 == 0) {
        for (auto &intersection : results) {
          std::cout << intersection.dist << std::endl;
          std::cout << intersection.path.x << "\t" << intersection.path.y << "\t" << intersection.path.z << std::endl;
        }
      }
      #endif
      results.clear();
  }
  check_sum = check_sum_v.sum();
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\tdone" << std::endl;
}

//
// Use a horizontal vectorization and data set
//
template <unsigned int kPlanes> void intersectVcHoriz(ray_data<Vector3<Scalar_v>>& ray,
                                                      plane_data<aligned::vector<vector_v>>& planes) {
  std::cout << "Running horiz. Vectoriztion test ...";
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  auto padding = planes.points.front().padding();

  // Access to raw data that will be loaded as scalar_v
  size_t n_bytes = planes.points.size() * (planes.points.front().n_elemts() + padding);
  size_t offset  = Scalar_v::Size + padding;
  size_t n_vec = n_bytes / (3*offset);

  BOOST_REQUIRE(n_bytes % (3*offset) == 0);

  for (size_t nt = 0; nt < nTests; ++nt) {
      auto pl_normals_ptr = const_cast<const Scalar*>(planes.normals.front().data());
      auto pl_points_ptr  = const_cast<const Scalar*>(planes.points.front().data());

      for (size_t i = 0; i < n_vec; i++) {
        auto intersection = vc_intersect_horiz<Scalar_v, const Scalar*>(ray, pl_normals_ptr, pl_points_ptr, offset);

        check_sum_v     += intersection.dist;
        pl_normals_ptr += 3 * offset;
        pl_points_ptr  += 3 * offset;

        #ifdef DEBUG 
        if (nt % (nTests-1)/2 == 0) {
          std::cout << intersection.dist << std::endl;
          std::cout << intersection.path.x << "\t" << intersection.path.y << "\t"
                    << intersection.path.z << std::endl;
        }
        #endif
      }
  }
  check_sum = check_sum_v.sum();
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\t\t\t\tdone" << std::endl;
}

//
// Use a horizontal vectorization and data set, save the results in a container
//
template <unsigned int kPlanes> void intersectVcHoriz_res(ray_data<Vector3<Scalar_v>>& ray,
                                                          plane_data<aligned::vector<vector_v>>& planes) {
  std::cout << "Running horiz. Vectoriztion with vector container test ...";
  Scalar   check_sum   = 0.0;
  Scalar_v check_sum_v = 0.0;

  auto padding = planes.points.front().padding();

  aligned::vector<intersection<Scalar_v, Vector3<Scalar_v>>> results;
  results.reserve(planes.points.size() * planes.points.front().n_elemts()/Scalar_v::Size);

  for (size_t nt = 0; nt < nTests; ++nt) {
      vc_intersect_horiz<Scalar_v, vector_v::obj_type>(ray, planes, results, padding);
      for (auto &intersection : results) check_sum_v += intersection.dist;
      #ifdef DEBUG 
      if (nt % (nTests-1)/2 == 0) {
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
  BOOST_CHECK(check_sum == -12*static_cast<int>(nTests));
  std::cout << "\tdone" << std::endl;
}


//----------------------------------------------------Run Boost nTests

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


