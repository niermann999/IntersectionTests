 /**
 * author: asalzburger@gmail.com
 **/

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include <Vc/Vc>

#include <array>
#include <iostream>
#include <limits>

namespace vecint {

namespace Test {


BOOST_AUTO_TEST_SUITE(VectIntersect)

unsigned int tests = 10000000;

using Scalar = float;

// Detray style intersection structure to be filled
template <typename scalar_t, typename vector_t>
struct intersection {
  scalar_t dist = std::numeric_limits<scalar_t>::infinity();
  vector_t path;
};

//=================================================================================================
// Eigen
//=================================================================================================

using Vector3_eig = Eigen::Matrix<Scalar, 3, 1>;
using Vector4_eig = Eigen::Vector4f;

//----------------------------------------------------Define Intersectors

auto intersect4D(Vector4_eig rayVector,
                 Vector4_eig rayPoint,
                 Vector4_eig planeNormal,
                 Vector4_eig planePoint) {
  float denom = rayVector.dot(planeNormal);
  float nom = (rayPoint - planePoint).dot(planeNormal);

  intersection<float, Vector4_eig> results = {.dist = nom / denom, 
                                              .path = Vector4_eig(rayPoint - nom*rayVector)};
  return results;
}

template <typename scalar_t, unsigned int kDIM>
auto intersect(Eigen::Matrix<scalar_t, kDIM, 3> rayVector,
               Eigen::Matrix<scalar_t, kDIM, 3> rayPoint,
               Eigen::Matrix<scalar_t, kDIM, 3> planeNormal,
               Eigen::Matrix<scalar_t, kDIM, 3> planePoint) {
  using Matrix_t = Eigen::Matrix<scalar_t, kDIM, 3>;
  using Array_t  = Eigen::Array<scalar_t, kDIM, 1>;
  
  Matrix_t tmpM_1 (std::move((rayPoint - planePoint).cwiseProduct(planeNormal)));
  Matrix_t tmpM_2 (std::move(rayVector.cwiseProduct(planeNormal)));
  Array_t coeffs ((tmpM_1.col(0) + tmpM_1.col(1) + tmpM_1.col(3)).array() 
                   / (tmpM_2.col(0) + tmpM_2.col(1) + tmpM_2.col(3)).array());

  // Broadcast coefficients onto ray-vector matrix
  rayVector.col(0).array() *= coeffs;
  rayVector.col(1).array() *= coeffs;
  rayVector.col(2).array() *= coeffs;
  
  intersection<Array_t, Matrix_t> results = {.dist = coeffs, 
                                             .path = rayPoint - rayVector};
  return results;
}

//----------------------------------------------------Fill Data Types

// Same starting position
Vector3_eig rv = Vector3_eig(0.0, -1.0, -1.0);
Vector3_eig rp = Vector3_eig(0.0, 0.0, 10.0);

// For the moment same normal vectors
Vector3_eig pn = Vector3_eig(0.0, 0.0, 1.0);

// 8 planes
Vector3_eig pp0_eg = Vector3_eig(0.0, 0.0, 5.0);
Vector3_eig pp1_eg = Vector3_eig(0.0, 0.0, 6.0);
Vector3_eig pp2_eg = Vector3_eig(0.0, 0.0, 7.0);
Vector3_eig pp3_eg = Vector3_eig(0.0, 0.0, 8.0);
Vector3_eig pp4_eg = Vector3_eig(0.0, 0.0, 9.0);
Vector3_eig pp5_eg = Vector3_eig(0.0, 0.0, 10.0);
Vector3_eig pp6_eg = Vector3_eig(0.0, 0.0, 11.0);
Vector3_eig pp7_eg = Vector3_eig(0.0, 0.0, 12.0);
Vector3_eig pp8_eg = Vector3_eig(0.0, 0.0, 13.0);

std::vector<Vector3_eig, Eigen::aligned_allocator<Vector3_eig> > pps_eg = {pp0_eg, pp1_eg, pp2_eg, pp3_eg, 
                                                                           pp4_eg, pp5_eg, pp6_eg, pp7_eg};

// Fixed size vectorizable
// Same starting position
Vector4_eig rv_4D = Vector4_eig(0.0, -1.0, -1.0, 0.0);
Vector4_eig rp_4D = Vector4_eig(0.0, 0.0, 10.0, 0.0);

// For the moment same normal vectors
Vector4_eig pn_4D = Vector4_eig(0.0, 0.0, 1.0, 0.0);

// 8 planes
Vector4_eig pp0_eg_4D = Vector4_eig(0.0, 0.0, 5.0, 0.0);
Vector4_eig pp1_eg_4D = Vector4_eig(0.0, 0.0, 6.0, 0.0);
Vector4_eig pp2_eg_4D = Vector4_eig(0.0, 0.0, 7.0, 0.0);
Vector4_eig pp3_eg_4D = Vector4_eig(0.0, 0.0, 8.0, 0.0);
Vector4_eig pp4_eg_4D = Vector4_eig(0.0, 0.0, 9.0, 0.0);
Vector4_eig pp5_eg_4D = Vector4_eig(0.0, 0.0, 10.0, 0.0);
Vector4_eig pp6_eg_4D = Vector4_eig(0.0, 0.0, 11.0, 0.0);
Vector4_eig pp7_eg_4D = Vector4_eig(0.0, 0.0, 12.0, 0.0);
Vector4_eig pp8_eg_4D = Vector4_eig(0.0, 0.0, 13.0, 0.0);

std::vector<Vector4_eig, Eigen::aligned_allocator<Vector4_eig> > pps_eg_4D = {pp0_eg_4D, pp1_eg_4D, pp2_eg_4D, 
                                                                              pp3_eg_4D, pp4_eg_4D, pp5_eg_4D, 
                                                                              pp6_eg_4D, pp7_eg_4D};

//----------------------------------------------------Run Tests

template <unsigned int kPlanes> void intersectSingle() {
  for (unsigned int nt = 0; nt < tests; ++nt) {
    for (unsigned int ip = 0; ip < kPlanes; ++ip) {
      auto intr = intersect<Scalar, 1>(rv, rp, pn, pps_eg[ip]);
      if (nt % 1000000 == 0) std::cout << intr.path << std::endl;
    }
  }
}

template <unsigned int kPlanes> void intersectSingle4D() {
  for (unsigned int nt = 0; nt < tests; ++nt) {
    for (unsigned int ip = 0; ip < kPlanes; ++ip) {
      auto intr = intersect4D(rv_4D, rp_4D, pn_4D, pps_eg_4D[ip]);
      if (nt % 1000000 == 0) std::cout << intr.path << std::endl;
    }
  }
}

template <unsigned int kPlanes> void intersectMultiple() {

  //constexpr unsigned int kFullDim = kPlanes * 3;

  using VectorM = Eigen::Matrix<Scalar, kPlanes, 3>;

  VectorM rvM;
  VectorM rpM;
  VectorM pnM;
  VectorM ppM;

  for (unsigned int ip = 0; ip < kPlanes; ++ip) {
    rvM.row(ip) << rv;
    rpM.row(ip) << rp;
    pnM.row(ip) << pn;
    ppM.row(ip) << pps_eg[ip];
  }

  std::vector<float> results;
  results.reserve(pps_eg.size());
  for (unsigned int nt = 0; nt < tests; ++nt) {
    auto intr = intersect<Scalar, kPlanes>(rvM, rpM, pnM, ppM);
    if (nt % 1000000 == 0) std::cout << intr.path(0) << std::endl;
  }
}

//=================================================================================================
// Vc
//=================================================================================================
using Scalar_v  = Vc::float_v;
using Vector3_v = Vc::float_v;
using Index_v   = Vector3_v::IndexType;

// Allow Vc to get alignment right
template <typename T, typename Allocator = Vc::Allocator<T>>
// Add subscript operator to allow for gather operations
using vector_aligned = Vc::Common::AdaptSubscriptOperator<std::vector<T, Allocator>>;

template<unsigned int kDIM>
using mem_t = Vc::Memory<Vector3_v, kDIM>;

// AoS
template<typename data_t>
struct Vector3
{
  data_t x, y, z;
};

// AoS for vertical vectorization
// Keep the geometrical vectors as Vc vectors (vertical vect.)
// Keep the plane points extra to make the struct alignment easier
template<typename data_t>
struct Vector3_vert
{
  data_t rayVector, rayPoint, planeNormal;
  //mem_t<kDIM> planePoints;
  //std::vector<data_t, Vc::Allocator<data_t>> planePoints;
};

//----------------------------------------------------Define Intersectors

template<typename vector_t, typename scalar_t, unsigned int kDIM>
auto vc_intersect_vert(Vector3_vert<vector_t> &data,
                       mem_t<kDIM> planePoints) {

  vector_t denom_v (data.rayVector * data.planeNormal);

  // Vector iterate
  auto iterations = planePoints.vectorsCount();
  std::vector<intersection<scalar_t, vector_t>> intersections;
  intersections.reserve(iterations);
  for (int i = 0; i < iterations; i++) {
    vector_t nom_v ((data.rayPoint - planePoints.vector(i)) * data.planeNormal);
    scalar_t coeff = (nom_v.sum() / denom_v.sum());
    intersections[i] = {.dist = coeff, 
                        .path = data.rayPoint - coeff*data.rayVector};
  }
  return intersections;
}


/*template<typename vector_t, typename scalar_t, unsigned int kDIM = 3>
auto vc_intersect_horiz(Vector3<vector_t> &rayVector,
                        Vector3<vector_t> &rayPoint,
                        Vector3<vector_t> &planeNormal,
                        std::vector<Vector3<scalar_t>> &pps_struct) {

  vector_t denom_x (rayVector.x * planeNormal.x);
  vector_t denom_y (rayVector.y * planeNormal.y);
  vector_t denom_z (rayVector.z * planeNormal.z);

  vector_t denoms (denom_x + denom_y + denom_z);

  vector_aligned<intersection<vector_t, Vector3<vector_t>>> intersections;
  intersections.reserve(pps_struct.size());
  int j = 0;
  // Vector iterate
  for (Index_v i(Vc::IndexesFromZero); (i < Index_v(pps_struct.size())).isFull(); i += Index_v(Vector3_v::Size)) {

    // Only possible without polymorphism etc.!
    Vector3_v x((vector_aligned<Vector3<Scalar>>::pointer)pps_struct.data(), &Vector3<scalar_t>::x, i);
    Vector3_v y((vector_aligned<Vector3<Scalar>>::pointer)pps_struct.data(), &Vector3<scalar_t>::y, i);
    Vector3_v z((vector_aligned<Vector3<Scalar>>::pointer)pps_struct.data(), &Vector3<scalar_t>::z, i);

    vector_t nom_x ((rayPoint.x - x) * planeNormal.x);
    vector_t nom_y ((rayPoint.y - y) * planeNormal.y);
    vector_t nom_z ((rayPoint.z - z) * planeNormal.z);

    vector_t coeffs ((nom_x + nom_y + nom_z)/denoms);
    Vector3<vector_t> path = {.x = (rayPoint.x - vector_t(coeffs[0])*rayVector.x),
                              .y = (rayPoint.y - vector_t(coeffs[1])*rayVector.y),
                              .z = (rayPoint.z - vector_t(coeffs[2])*rayVector.z)};
    intersections[j++] = {.dist = coeffs, .path = path};
  }
  return intersections;
}*/

template<typename vector_t, unsigned int kDIM = 3>
auto vc_intersect_horiz(Vector3<vector_t> &rayVector,
                        Vector3<vector_t> &rayPoint,
                        Vector3<vector_t> &planeNormal,
                        Vector3<mem_t<kDIM> > &planePoints) {
  
  auto iters = planePoints.x.vectorsCount();
  vector_aligned<intersection<vector_t, Vector3<vector_t>>> intersections;
  intersections.reserve(iters);

  for (int i = 0; i < iters; i++) {
    vector_t denom_x (rayVector.x * planeNormal.x);
    vector_t denom_y (rayVector.y * planeNormal.y);
    vector_t denom_z (rayVector.z * planeNormal.z);

    vector_t denoms (denom_x + denom_y + denom_z);

    vector_t nom_x ((rayPoint.x - planePoints.x.vector(i, Vc::Streaming)) * planeNormal.x);
    vector_t nom_y ((rayPoint.y - planePoints.y.vector(i, Vc::Streaming)) * planeNormal.y);
    vector_t nom_z ((rayPoint.z - planePoints.z.vector(i, Vc::Streaming)) * planeNormal.z);

    vector_t coeffs ((nom_x + nom_y + nom_z)/denoms);

    Vector3<vector_t> path = {.x = (rayPoint.x - vector_t(coeffs[0])*rayVector.x),
                              .y = (rayPoint.y - vector_t(coeffs[1])*rayVector.y),
                              .z = (rayPoint.z - vector_t(coeffs[2])*rayVector.z)};
    intersections[i] = {.dist = coeffs, .path = path};
  }
  return intersections;
}

template<typename vector_t, unsigned int kDIM = 3>
auto vc_intersect_horiz(Vector3<vector_t> &rayVector,
                        Vector3<vector_t> &rayPoint,
                        Vector3<vector_t> &planeNormal,
                        mem_t<kDIM> &planePoints) {

  /*vector_t denoms (rayVector.x*planeNormal.x + rayVector.y*planeNormal.y + rayVector.z*planeNormal.z);

  auto iters = planePoints.vectorsCount()/3;
  auto mem_itr = planePoints.begin(Vc::Streaming);
  std::vector<intersection<vector_t, Vector3<vector_t>>> intersections;
  intersections.reserve(iters);
  for (int i = 0; i < iters; i++) {
    vector_t coeffs ((rayPoint.x - *std::next(mem_itr)) * planeNormal.x
                    +(rayPoint.y - *std::next(mem_itr)) * planeNormal.y
                    +(rayPoint.z - *std::next(mem_itr)) * planeNormal.z);
    coeffs /= denoms;

    Vector3<vector_t> path = {.x = rayPoint.x - rayVector.x*coeffs, 
                              .y = rayPoint.y - rayVector.y*coeffs, 
                              .z = rayPoint.z - rayVector.z*coeffs};
    intersections[i] = {.dist = std::move(coeffs), .path = std::move(path)};
  }
  return std::move(intersections);*/
  intersection<vector_t, Vector3<vector_t>> inter = {.dist=rayVector.x, .path=rayPoint};
  std::vector<intersection<vector_t, Vector3<vector_t>>> inters;
  inters.push_back(inter);
  //inters[0].dist += kDIM;
  return inters;
}

//----------------------------------------------------Fill Data Types

//**************************AoS

Vector3<Scalar> plain_normal {.x= 0.0, .y=0.0,  .z=1.0};
Vector3<Scalar> ray_vector   {.x= 0.0, .y=-1.0, .z=-1.0};
Vector3<Scalar> ray_point    {.x= 0.0, .y=0.0,  .z=10.0};

Vector3<Scalar> pp0 {.x=0.0, .y=0.0, .z=5.0};
Vector3<Scalar> pp1 {.x=0.0, .y=0.0, .z=6.0};
Vector3<Scalar> pp2 {.x=0.0, .y=0.0, .z=7.0};
Vector3<Scalar> pp3 {.x=0.0, .y=0.0, .z=8.0};
Vector3<Scalar> pp4 {.x=0.0, .y=0.0, .z=9.0};
Vector3<Scalar> pp5 {.x=0.0, .y=0.0, .z=10.0};
Vector3<Scalar> pp6 {.x=0.0, .y=0.0, .z=11.0};
Vector3<Scalar> pp7 {.x=0.0, .y=0.0, .z=12.0};
Vector3<Scalar> pp8 {.x=0.0, .y=0.0, .z=13.0};

std::vector<Vector3<Scalar>> pps_struct = {pp0, pp1, pp2, pp3, pp4, pp5, pp6, pp7};
// Translate the Eigen data structure
mem_t<24> planePoints;

void fill_vert_Vec() {
  // vector based access:
  for (size_t i = 0; i < planePoints.vectorsCount(); i++) {
    auto pp = Vector3_v((pps_eg.data())[Vector3_v::Size * i]);
    planePoints.vector(i) = pp;
  }
}

Vector3_vert<Vector3_v> intersection_data = {.rayVector   = Vector3_v(rv.array().data()),
                                             .rayPoint    = Vector3_v(rp.array().data()), 
                                             .planeNormal = Vector3_v(pn.array().data())};

//**************************SoA
// As structs of vectors
Vector3<Scalar_v> pp_v;
Vector3<Scalar_v> pn_v {.x= Vector3_v(0.0), .y=Vector3_v(0.0),  .z=Vector3_v(1.0)};
Vector3<Scalar_v> rv_v {.x= Vector3_v(0.0), .y=Vector3_v(-1.0), .z=Vector3_v(-1.0)};
Vector3<Scalar_v> rp_v {.x= Vector3_v(0.0), .y=Vector3_v(0.0),  .z=Vector3_v(10.0)};

// As memory block
Vector3<mem_t<8>> pp_mem;

void fill_SoA () {
  mem_t<8> pp_mem_x;
  mem_t<8> pp_mem_y;
  mem_t<8> pp_mem_z;
  Vector3_v z1 = Vector3_v(Vc::IndexesFromZero);
  Vector3_v z2 = Vector3_v(Vc::IndexesFromZero);
  z1 += Vector3_v(5.0);
  z2 += Vector3_v(9.0);

  pp_mem_x.vector(0) = Vector3_v(0.0);
  pp_mem_x.vector(1) = Vector3_v(0.0);
  pp_mem_y.vector(0) = Vector3_v(0.0);
  pp_mem_y.vector(1) = Vector3_v(0.0);
  pp_mem_z.vector(0) = z1;
  pp_mem_z.vector(1) = z2;

  pp_mem.x = pp_mem_x;
  pp_mem.y = pp_mem_y;
  pp_mem.z = pp_mem_z;
}

// Fill everything into single memory interleaved as vc vectors to optimize memory loacality
mem_t<24> pps_mem;

void fill_inteleaved_mem () {
  Vector3_v z1 = Vector3_v(Vc::IndexesFromZero);
  Vector3_v z2 = Vector3_v(Vc::IndexesFromZero);
  z1 += Vector3_v(5.0);
  z2 += Vector3_v(9.0);

  pps_mem.vector(0) = Vector3_v(0.0);
  pps_mem.vector(3) = Vector3_v(0.0);
  pps_mem.vector(1) = Vector3_v(0.0);
  pps_mem.vector(4) = Vector3_v(0.0);
  pps_mem.vector(2) = z1;
  pps_mem.vector(5) = z2;
}

//----------------------------------------------------Run Tests

void intersectVc_vert() {
    for (unsigned int nt = 0; nt < tests; ++nt) {
      auto intr = vc_intersect_vert<Vector3_v, Scalar, 24>(intersection_data, planePoints);
      if (nt % 1000000 == 0) std::cout << intr[0].dist << std::endl;
    }
}

/*void intersectVc_horiz_AoS() {
    for (unsigned int nt = 0; nt < tests; ++nt) {
      auto intr = vc_intersect_horiz<Vector3_v, Scalar>(rv_v, rp_v, pn_v, pps_struct);
      if (nt % 100000 == 0) std::cout << intr[0].dist << std::endl;
    }
}*/
  

void intersectVc_horiz_SoA() {
    for (unsigned int nt = 0; nt < tests; ++nt) {
      auto intr = vc_intersect_horiz<Vector3_v, 8>(rv_v, rp_v, pn_v, pp_mem);
      if (nt % 100000 == 0) std::cout << intr[0].dist << std::endl;
    }
}

void intersectVc_horiz_SoA_interleaved() {
    for (unsigned int nt = 0; nt < tests; ++nt) {
      auto intr = vc_intersect_horiz<Vector3_v, 24>(rv_v, rp_v, pn_v, pps_mem);
      //if (nt % 100000 == 0) std::cout << intr[0].dist << "\t" << intr[0].path.x << std::endl;
      if (nt % 100000 == 0) std::cout << intr[0].dist[0] << std::endl;
    }
}


//----------------------------------------------------Run Boost Tests

//BOOST_AUTO_TEST_CASE(SingleIntersection4) { intersectSingle<4>(); }

//BOOST_AUTO_TEST_CASE(MultiIntersection4) { intersectMultiple<4>(); }

//BOOST_AUTO_TEST_CASE(SingleIntersection6) { intersectSingle<6>(); }

//BOOST_AUTO_TEST_CASE(MultiIntersection6) { intersectMultiple<6>(); }

BOOST_AUTO_TEST_CASE(SingleIntersection8) { intersectSingle<8>(); }

//BOOST_AUTO_TEST_CASE(SingleIntersection4D8) { intersectSingle4D<8>(); }

//BOOST_AUTO_TEST_CASE(MultiIntersection8) { intersectMultiple<8>(); }

//BOOST_AUTO_TEST_CASE(SingleIntersection9) { intersectSingle<9>(); }

//BOOST_AUTO_TEST_CASE(MultiIntersection9) { intersectMultiple<9>(); }



BOOST_AUTO_TEST_CASE(fillSoA) { fill_SoA(); }

BOOST_AUTO_TEST_CASE(fillSoA_inter) { fill_inteleaved_mem(); }

BOOST_AUTO_TEST_CASE(fill_vert) { fill_vert_Vec(); }

BOOST_AUTO_TEST_CASE(VcIntersect_Vert) { intersectVc_vert(); }

//BOOST_AUTO_TEST_CASE(VcIntersect_AoS) { intersectVc_horiz_AoS(); }

//BOOST_AUTO_TEST_CASE(VcIntersect_SoA) { intersectVc_horiz_SoA(); }

BOOST_AUTO_TEST_CASE(VcIntersect_SoA_inter) { intersectVc_horiz_SoA_interleaved(); }



BOOST_AUTO_TEST_SUITE_END()

} // namespace Test

} // namespace vecint


