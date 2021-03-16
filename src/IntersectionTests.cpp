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

unsigned int tests = 1000000;

using Scalar = float;

//=================================================================================================
// Eigen
//=================================================================================================

using Vector3_eig = Eigen::Matrix<Scalar, 3, 1>;
using Vector4_eig = Eigen::Vector4f;

//----------------------------------------------------Define Intersectors

auto intersect4D(Vector4_eig rayVector,
                 Vector4_eig rayPoint,
                 Vector4_eig planeNormal,
                 Vector4_eig planePoint,
                 std::vector<float> &results) {
  float denom = rayVector.dot(planeNormal);
  float nom = (rayPoint - planePoint).dot(planeNormal);
  results.push_back(nom / denom);
}

template <typename scalar_t, unsigned int kDIM>
auto intersect(Eigen::Matrix<scalar_t, kDIM, 3> rayVector,
               Eigen::Matrix<scalar_t, kDIM, 3> rayPoint,
               Eigen::Matrix<scalar_t, kDIM, 3> planeNormal,
               Eigen::Matrix<scalar_t, kDIM, 3> planePoint,
               std::vector<float> &results) {
  Eigen::Matrix<scalar_t, kDIM, 3> tmpM_1 (std::move((rayPoint - planePoint).cwiseProduct(planeNormal)));
  Eigen::Matrix<scalar_t, kDIM, 3> tmpM_2 (std::move(rayVector.cwiseProduct(planeNormal)));
  Eigen::Array<scalar_t, kDIM, 1> coeffs ((tmpM_1.col(0) + tmpM_1.col(1) + tmpM_1.col(3)).array() 
                                           / (tmpM_2.col(0) + tmpM_2.col(1) + tmpM_2.col(3)).array());
  
  /*// Broadcast coefficients onto ray-vector matrix
  rayVector.col(0).array() *= coeffs;
  rayVector.col(1).array() *= coeffs;
  rayVector.col(2).array() *= coeffs;
  
  return rayPoint - rayVector;*/
  //return coeffs;
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
    std::vector<float> results;
    results.reserve(kPlanes);
    for (unsigned int ip = 0; ip < kPlanes; ++ip) {
      intersect<Scalar, 1>(rv, rp, pn, pps_eg[ip], results);
      if (nt % 100000 == 0) std::cout << results[0] << std::endl;
    }
  }
}

template <unsigned int kPlanes> void intersectSingle4D() {
  std::vector<float> results;
  results.reserve(pps_eg.size());
  for (unsigned int nt = 0; nt < tests; ++nt) {
    for (unsigned int ip = 0; ip < kPlanes; ++ip) {
      intersect4D(rv_4D, rp_4D, pn_4D, pps_eg_4D[ip], results);
      if (nt % 100000 == 0) std::cout << results[0] << std::endl;
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
    intersect<Scalar, kPlanes>(rvM, rpM, pnM, ppM, results);
    if (nt % 100000 == 0) std::cout << results[0] << std::endl;
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
template<typename data_t, unsigned int kDIM>
struct Vector3_vert
{
  data_t rayVector, rayPoint, planeNormal;
  //mem_t<kDIM> planePoints;
  //std::vector<data_t, Vc::Allocator<data_t>> planePoints;
};

//----------------------------------------------------Define Intersectors

template<typename scalar_t, typename vector_t, unsigned int kDIM>
auto vc_intersect_vert(Vector3_vert<vector_t, kDIM> &data,
                       mem_t<kDIM> planePoints,
                       vector_aligned<scalar_t> &results) {

  vector_t denom_v (data.rayVector * data.planeNormal);

  // Vector iterate
  for (int i = 0; i < planePoints.vectorsCount(); i++) {
    vector_t nom_v ((data.rayPoint - planePoints.vector(i)) * data.planeNormal);
    results[i] = (nom_v.sum() / denom_v.sum());
  }
}


template<typename data_t, unsigned int kDIM = 3>
auto vc_intersect_horiz(Vector3<data_t> &rayVector,
                        Vector3<data_t> &rayPoint,
                        Vector3<data_t> &planeNormal,
                        std::vector<Vector3<Scalar>> &pps_struct,
                        vector_aligned<data_t> &results) {

  data_t denom_x (rayVector.x * planeNormal.x);
  data_t denom_y (rayVector.y * planeNormal.y);
  data_t denom_z (rayVector.z * planeNormal.z);

  data_t denoms (denom_x + denom_y + denom_z);
  int j = 0;
  // Vector iterate
  for (Index_v i(Vc::IndexesFromZero); (i < Index_v(pps_struct.size())).isFull(); i += Index_v(Vector3_v::Size)) {

    // Only possible without polymorphism etc.!
    Vector3_v x((vector_aligned<Vector3<Scalar>>::pointer)pps_struct.data(), &Vector3<Scalar>::x, i);
    Vector3_v y((vector_aligned<Vector3<Scalar>>::pointer)pps_struct.data(), &Vector3<Scalar>::y, i);
    Vector3_v z((vector_aligned<Vector3<Scalar>>::pointer)pps_struct.data(), &Vector3<Scalar>::z, i);

    data_t nom_x ((rayPoint.x - x) * planeNormal.x);
    data_t nom_y ((rayPoint.y - y) * planeNormal.y);
    data_t nom_z ((rayPoint.z - z) * planeNormal.z);

    data_t noms (nom_x + nom_y + nom_z);

    results[j++] = noms / denoms;
  }
}

template<typename data_t, unsigned int kDIM = 3>
auto vc_intersect_horiz(Vector3<data_t> &rayVector,
                        Vector3<data_t> &rayPoint,
                        Vector3<data_t> &planeNormal,
                        Vector3<mem_t<kDIM> > &planePoints,
                        vector_aligned<data_t> &results) {

  data_t denom_x (rayVector.x * planeNormal.x);
  data_t denom_y (rayVector.y * planeNormal.y);
  data_t denom_z (rayVector.z * planeNormal.z);
  
  for (int i = 0; i < 2; i++) {
    data_t nom_x ((rayPoint.x - planePoints.x.vector(i)) * planeNormal.x);
    data_t nom_y ((rayPoint.y - planePoints.y.vector(i)) * planeNormal.y);
    data_t nom_z ((rayPoint.z - planePoints.z.vector(i)) * planeNormal.z);

    data_t noms (nom_x + nom_y + nom_z);
    data_t denoms (denom_x + denom_y + denom_z);

    results[i] = noms / denoms;
  }
}

template<typename data_t, unsigned int kDIM = 3>
auto vc_intersect_horiz(Vector3<data_t> &rayVector,
                        Vector3<data_t> &rayPoint,
                        Vector3<data_t> &planeNormal,
                        mem_t<kDIM> &planePoints,
                        vector_aligned<data_t> &results) {

  data_t denom_x (rayVector.x * planeNormal.x);
  data_t denom_y (rayVector.y * planeNormal.y);
  data_t denom_z (rayVector.z * planeNormal.z);

  data_t denoms (denom_x + denom_y + denom_z);

  int j = 0;
  for (int i = 0; i < 2; i++) {
    data_t nom_x ((rayPoint.x - planePoints.vector(j++)) * planeNormal.x);
    data_t nom_y ((rayPoint.y - planePoints.vector(j++)) * planeNormal.y);
    data_t nom_z ((rayPoint.z - planePoints.vector(j++)) * planeNormal.z);

    data_t noms (nom_x + nom_y + nom_z);

    results[i] = noms / denoms;
  }
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

Vector3_vert<Vector3_v, 24> intersection_data = {.rayVector   = Vector3_v(rv.array().data()),
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

    vector_aligned<Scalar> results;
    results.reserve(pps_eg.size());
    for (unsigned int nt = 0; nt < tests; ++nt) {
      vc_intersect_vert<Scalar, Vector3_v, 24>(intersection_data, planePoints, results);
      if (nt % 100000 == 0) std::cout << results[0] << "\n" << results[1] << std::endl;
    }
}

void intersectVc_horiz_AoS() {

    vector_aligned<Vector3_v> results;
    results.reserve(2*Vector3_v::Size);
    for (unsigned int nt = 0; nt < tests; ++nt) {
      vc_intersect_horiz<Vector3_v>(rv_v, rp_v, pn_v, pps_struct, results);
      if (nt % 100000 == 0) std::cout << results[0] << "\n" << results[1] << std::endl;
    }
}
  

void intersectVc_horiz_SoA() {

    vector_aligned<Vector3_v> results;
    results.reserve(2*Vector3_v::Size);
    for (unsigned int nt = 0; nt < tests; ++nt) {
      vc_intersect_horiz<Vector3_v, 8>(rv_v, rp_v, pn_v, pp_mem, results);
      if (nt % 100000 == 0) std::cout << results[0] << "\n" << results[1] << std::endl;
    }
}

void intersectVc_horiz_SoA_interleaved() {
  
    vector_aligned<Vector3_v> results;
    results.reserve(2*Vector3_v::Size);
    for (unsigned int nt = 0; nt < tests; ++nt) {
      vc_intersect_horiz<Vector3_v, 24>(rv_v, rp_v, pn_v, pps_mem, results);
      if (nt % 100000 == 0) std::cout << results[0] << "\n" << results[1] << std::endl;
    }
}


//----------------------------------------------------Run Boost Tests

//BOOST_AUTO_TEST_CASE(SingleIntersection4) { intersectSingle<4>(); }

//BOOST_AUTO_TEST_CASE(MultiIntersection4) { intersectMultiple<4>(); }

//BOOST_AUTO_TEST_CASE(SingleIntersection6) { intersectSingle<6>(); }

//BOOST_AUTO_TEST_CASE(MultiIntersection6) { intersectMultiple<6>(); }

BOOST_AUTO_TEST_CASE(SingleIntersection8) { intersectSingle<8>(); }

BOOST_AUTO_TEST_CASE(SingleIntersection4D8) { intersectSingle4D<8>(); }

BOOST_AUTO_TEST_CASE(MultiIntersection8) { intersectMultiple<8>(); }

//BOOST_AUTO_TEST_CASE(SingleIntersection9) { intersectSingle<9>(); }

//BOOST_AUTO_TEST_CASE(MultiIntersection9) { intersectMultiple<9>(); }



BOOST_AUTO_TEST_CASE(fillSoA) { fill_SoA(); }

BOOST_AUTO_TEST_CASE(fillSoA_inter) { fill_inteleaved_mem(); }

BOOST_AUTO_TEST_CASE(fill_vert) { fill_vert_Vec(); }

BOOST_AUTO_TEST_CASE(VcIntersect_Vert) { intersectVc_vert(); }

BOOST_AUTO_TEST_CASE(VcIntersect_AoS) { intersectVc_horiz_AoS(); }

BOOST_AUTO_TEST_CASE(VcIntersect_SoA) { intersectVc_horiz_SoA(); }

BOOST_AUTO_TEST_CASE(VcIntersect_SoA_inter) { intersectVc_horiz_SoA_interleaved(); }



BOOST_AUTO_TEST_SUITE_END()

} // namespace Test

} // namespace vecint
