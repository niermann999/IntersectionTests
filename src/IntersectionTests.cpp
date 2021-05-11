 /**
 * author: joana.niermann@cern.ch
 **/
#include "types.hpp"
#include "intersectors.hpp"

#include <atomic>
//#include <barrier> //c++ 20
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

using vector_s = data_trait<Vector4_s>::type;
using vector_v = data_trait<Vector3_v>::type;

//----------------------------------------------------Fill Data Types

// Produce test data
struct data_setup {

  vector_v ray_dir_hor, ray_point_hor;

  vector_s        pl_normals, pl_points;
  aligned::vector<Vector3<Scalar> > pl_points_struct, pl_normals_struct;
  vector_v        pl_normals_hor, pl_points_hor;

  ray_data<vector_s> ray;
  aligned::vector<plane_data<vector_s> > planes;

  ray_data<vector_v> ray_struct;
  plane_data<aligned::vector<Vector3<Scalar> > > planes_struct;

  ray_data<vector_v> ray_hor;
  aligned::vector<plane_data<vector_v> > planes_hor;

  data_setup () {
    // Ray to be intersected with different plains
    Scalar data[4] = {0.0, -1.0, -1.0, 0.0};//, 0.0, 0.0, 0.0, 0.0, 0.0};
    vector_s ray_dir {.obj = vector_s::type(data)};
    data[2] = 10.0;
    data[1] = 0.0;
    vector_s ray_point {.obj = vector_s::type(data)};
    data[2] = 1.0;

    // For the moment same normal vectors
    vector_s pl_normal {.obj = vector_s::type(data)};

    // plane normals and points
    data[2] = 5.0;
    planes.reserve(nSurfaces);
    for (size_t i = 0; i < nSurfaces; i++) {
      vector_s pl_point {.obj = vector_s::type(data)};
      planes.push_back({.normals = pl_normal, .points = pl_point});
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
    pl_normals_struct = {plain_n, plain_n, plain_n, plain_n, plain_n, plain_n, plain_n, plain_n};
    planes_struct = {.normals = pl_normals_struct, .points  = pl_points_struct};

    // horizontal data (interleaved)
    vector_v ray_dir_hor, ray_point_hor;
    ray_dir_hor.obj   = {.x= Scalar_v(0.0), .y=Scalar_v(-1.0), .z=Scalar_v(-1.0)};
    ray_point_hor.obj = {.x= Scalar_v(0.0), .y=Scalar_v(0.0), .z=Scalar_v(10.0)};


    for (size_t i = 0; i < planes.size()/Scalar_v::Size; i++) {
      for (size_t j = 0; j < Scalar_v::Size; j++) {
        pl_points_hor.obj.x[j]  = planes[j+i*Scalar_v::Size].points()[0];
        pl_normals_hor.obj.x[j] = planes[j+i*Scalar_v::Size].normals()[0];
        pl_points_hor.obj.y[j]  = planes[j+i*Scalar_v::Size].points()[1];
        pl_normals_hor.obj.y[j] = planes[j+i*Scalar_v::Size].normals()[1];
        pl_points_hor.obj.z[j]  = planes[j+i*Scalar_v::Size].points()[2];
        pl_normals_hor.obj.z[j] = planes[j+i*Scalar_v::Size].normals()[2];
      }
      planes_hor.push_back({.normals = pl_normals_hor, .points = pl_points_hor});
    }

    // vertical data containers
    ray    = {.direction = ray_dir, 
              .point     = ray_point};

    // AoS container
    ray_struct = {.direction = ray_dir_hor, 
                  .point     = ray_point_hor};

    // horizontal ray container
    ray_hor = {.direction = ray_dir_hor,    
               .point     = ray_point_hor};
  }
};

//----------------------------------------------------Run Tests

//----------------------------------------------------//
// Eigen                                              //
//----------------------------------------------------//
template <unsigned int kPlanes> void intersectEigen4D(ray_data<vector_s>& ray,
                                                      aligned::vector<plane_data<vector_s> >& planes) {
  using scalar_t = typename vector_s::scalar_type;

  std::cout << "Running Eigen 4D test ...";
  scalar_t check_sum = 0.0;

  for (auto &plane : planes) {
    auto intersection = eig_intersect_4D(ray, plane);
    check_sum += intersection.dist;
    #ifdef DEBUG
    std::cout << "\n" << intersection.dist << std::endl;
    for (int i = 0; i < 4; i++) {
      std::cout << intersection.path[i] << ", ";
    }
    #endif
  }

  // TODO watch the floating point comparision!!
  BOOST_CHECK(check_sum == 12);
  std::cout << "\t\t\t\t\tdone" << std::endl;
}

template <unsigned int kPlanes> void intersectEigen4D_res(ray_data<vector_s>& ray,
                                                          aligned::vector<plane_data<vector_s> >& planes) {
  using scalar_t = typename vector_s::scalar_type;
  using vector_t = typename vector_s::type;

  std::cout << "Running Eigen with vector container 4D test ...";
  scalar_t check_sum = 0.0;

  aligned::vector<intersection<scalar_t, vector_t> > results;
  results.reserve(planes.size());

  eig_intersect_4D<vector_s>(ray, planes, results);
  for (auto &intersection : results) check_sum += intersection.dist;
  #ifdef DEBUG
  for (auto &intersection : results) {
    std::cout << "\n" << intersection.dist << std::endl;
    for (int i = 0; i < 4; i++) {
      std::cout << intersection.path[i] << ", ";
    }
  }
  #endif
  results.clear();

  BOOST_CHECK(check_sum == 12);
  std::cout << "\t\t\tdone" << std::endl;
}

//----------------------------------------------------//
// Vc                                                 //
//----------------------------------------------------//

//
// Use simdArray on unchanged data set
//
template <unsigned int kPlanes> void intersectVcVert(ray_data<vector_s>& ray,
                                                     aligned::vector<plane_data<vector_s> >& planes) {
  using scalar_t = typename vector_s::scalar_type;

  std::cout << "Running vert. Vectoriztion test ...";
  scalar_t check_sum = 0.0;

  for (auto &plane : planes) {
    auto intersection = vc_intersect_vert<vector_s>(ray, plane);
    check_sum += intersection.dist;
    #ifdef DEBUG
    std::cout << "\n" << intersection.dist << std::endl;
    for (int i = 0; i < 4; i++) {
      std::cout << intersection.path[i] << ", ";
    }
    #endif
  }
  BOOST_CHECK(check_sum == 12);
  std::cout << "\t\t\t\tdone" << std::endl;
}

//
// Use simdArray on unchanged data set and save the results in a container
//
template <unsigned int kPlanes> void intersectVcVert_res(ray_data<vector_s> &ray,
                                                         aligned::vector<plane_data<vector_s> > &planes) {
  using scalar_t = typename vector_s::scalar_type;

  std::cout << "Running vert. Vectoriztion with vector container test ...";
  scalar_t check_sum = 0.0;
  
  aligned::vector<intersection<scalar_t, Vc::SimdArray<scalar_t, 4> >> results;
  results.reserve(planes.size());

  vc_intersect_vert<vector_s>(ray, planes, results);
  for (auto &intersection : results) check_sum += intersection.dist;
  #ifdef DEBUG
  for (auto &intersection : results) {
    std::cout << "\n" << intersection.dist << std::endl;
    for (int i = 0; i < 4; i++) {
      std::cout << intersection.path[i] << ", ";
    }
  }
  #endif
  results.clear();
  
  BOOST_CHECK(check_sum == 12);
  std::cout << "\tdone" << std::endl;
}

//
// Use a gather on a structured data set then vectorize horizontaly
//
template <unsigned int kPlanes> void intersectVcHybrid(ray_data<vector_v> &ray,
                                                       plane_data<aligned::vector<Vector3<Scalar> >> &planes) {
  using scalar_t = typename vector_v::scalar_type;
  using scalar_v = typename vector_v::vec_type;

  std::cout << "Running hybrid Vectoriztion test ...";
  scalar_t   check_sum   = 0.0;
  scalar_v check_sum_v = 0.0;

  for (Index_v i(scalar_v::IndexesFromZero()); (i < Index_v(planes.points.size())).isFull(); i += Index_v(scalar_v::Size)) {
    vector_v pl_normal_strc, pl_point_strc;

    scalar_v pns_x = planes.normals[i][&Vector3<scalar_t>::x];
    scalar_v pns_y = planes.normals[i][&Vector3<scalar_t>::y];
    scalar_v pns_z = planes.normals[i][&Vector3<scalar_t>::z];
    pl_normal_strc.obj = {.x = pns_x, .y = pns_y, .z = pns_z};

    scalar_v pps_x = planes.points[i][&Vector3<scalar_t>::x];
    scalar_v pps_y = planes.points[i][&Vector3<scalar_t>::y];
    scalar_v pps_z = planes.points[i][&Vector3<scalar_t>::z];
    pl_point_strc.obj = {.x = pps_x, .y = pps_y, .z = pps_z};
    plane_data<vector_v> planes_strcts = {.normals = pl_normal_strc, .points = pl_point_strc};

    auto intersection = vc_intersect_hybrid<vector_v>(ray, planes_strcts);
    check_sum_v += intersection.dist;
    #ifdef DEBUG
      std::cout << intersection.dist << std::endl;
      std::cout << intersection.path.x << "\t" << intersection.path.y << "\t" << intersection.path.z << std::endl;
    #endif
  }
  check_sum = check_sum_v.sum();
  BOOST_CHECK(check_sum == 12);
  std::cout << "\t\t\t\tdone" << std::endl;
}

//
// Use a gather on a structured data set then vectorize horizontaly 
// and save the results in a container
//
/*template <unsigned int kPlanes> void intersectVcHybrid_res(ray_data<vector_v> &ray,
                                                           plane_data<aligned::vector<Vector3<Scalar> >> &planes) {
  using scalar_t = typename vector_v::scalar_type;
  using scalar_v = typename vector_v::vec_type;

  std::cout << "Running hybrid Vectoriztion with vector container test ...";
  scalar_t   check_sum   = 0.0;
  scalar_v check_sum_v = 0.0;

  aligned::vector<intersection<scalar_v, Vector3<Scalar_v> >> results;
  results.reserve(planes.normals.size());

  vc_intersect_hybrid<vector_v>(ray, planes, results);
  for (auto &intersection : results) check_sum_v += intersection.dist;
  #ifdef DEBUG
  for (auto &intersection : results) {
    std::cout << intersection.dist << std::endl;
    std::cout << intersection.path.x << "\t" << intersection.path.y << "\t" << intersection.path.z << std::endl;
  }
  #endif
  results.clear();

  check_sum = check_sum_v.sum();
  BOOST_CHECK(check_sum == 12);
  std::cout << "\tdone" << std::endl;
}*/

//
// Use a horizontal vectorization and data set
//
template <unsigned int kPlanes> void intersectVcHoriz(ray_data<vector_v> &ray,
                                                      aligned::vector<plane_data<vector_v> > &planes) {
  using scalar_t = typename vector_v::scalar_type;
  using scalar_v = typename vector_v::vec_type;

  std::cout << "Running horiz. Vectoriztion test ...";
  scalar_t check_sum   = 0.0;
  scalar_v check_sum_v = 0.0;

  for (auto &plane : planes) {
    auto intersection = vc_intersect_horiz<vector_v>(ray, plane);

    check_sum_v += intersection.dist;

    #ifdef DEBUG
    std::cout << intersection.dist << std::endl;
    std::cout << intersection.path.x << "\t" << intersection.path.y << "\t"
              << intersection.path.z << std::endl;
    #endif
  }
  check_sum = check_sum_v.sum();
  BOOST_CHECK(check_sum == 12);
  std::cout << "\t\t\t\tdone" << std::endl;
}

//
// Use a horizontal vectorization and data set, save the results in a container
//
template <unsigned int kPlanes> void intersectVcHoriz_res(ray_data<vector_v> &ray,
                                                          aligned::vector<plane_data<vector_v> > &planes) {
  using scalar_t = typename vector_v::scalar_type;
  using scalar_v = typename vector_v::vec_type;
  using vector_t = typename vector_v::type;

  std::cout << "Running horiz. Vectoriztion with vector container test ...";
  scalar_t check_sum   = 0.0;
  scalar_v check_sum_v = 0.0;

  aligned::vector<intersection<scalar_v, vector_t> > results;
  results.reserve(planes.size());

  vc_intersect_horiz<vector_v>(ray, planes, results);

  for (auto &intersection : results) check_sum_v += intersection.dist;
  #ifdef DEBUG
  for (auto &intersection : results) {
    std::cout << intersection.dist << std::endl;
    std::cout << intersection.path.x << "\t" << intersection.path.y << "\t"
              << intersection.path.z << std::endl;
  }
  #endif
  results.clear();

  check_sum = check_sum_v.sum();
  //BOOST_CHECK(check_sum == 12);
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
                                                                          
/*BOOST_AUTO_TEST_CASE(IntersectVcHybrid_res) {intersectVcHybrid_res<nSurfaces>(ray_struct, 
                                                                              planes_struct);}*/


BOOST_AUTO_TEST_CASE(IntersectVcHoriz)    {intersectVcHoriz<nSurfaces>(ray_hor, 
                                                                       planes_hor);}
BOOST_AUTO_TEST_CASE(IntersectVcHorizres) {intersectVcHoriz_res<nSurfaces>(ray_hor, 
                                                                           planes_hor);}

BOOST_AUTO_TEST_SUITE_END()

} // namespace test

} //namespace vec_intr


