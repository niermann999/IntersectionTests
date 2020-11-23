 /**
 * author: asalzburger@gmail.com
 **/

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include <array>
#include <iostream>
#include <limits>

namespace vecint {

namespace Test {

template <typename scalar_t, unsigned int kDIM>
auto intersect(Eigen::Matrix<scalar_t, kDIM, 1> rayVector,
               Eigen::Matrix<scalar_t, kDIM, 1> rayPoint,
               Eigen::Matrix<scalar_t, kDIM, 1> planeNormal,
               Eigen::Matrix<scalar_t, kDIM, 1> planePoint) {
  /*Eigen::Matrix<scalar_t, kDIM, 1> diff = rayPoint - planePoint;
  scalar_t prod1 = diff.dot(planeNormal);
  scalar_t prod2 = rayVector.dot(planeNormal);
  scalar_t prod3 = prod1 / prod2;
  return rayPoint - rayVector * prod3;*/
  
  //typedef Eigen::Array <scalar_t, kDIM/3, 1> array_t;
  //typedef Eigen::Matrix<scalar_t, kDIM,   1> vector_t;
  //typedef Eigen::Matrix<scalar_t, 3, kDIM/3> matrix_t;  
  
  //Eigen::Map<matrix_t> map1 (((vector_t)(rayPoint - planePoint).cwiseProduct(planeNormal)).data());
  //Eigen::Map<matrix_t> map2 (((vector_t)rayVector.cwiseProduct(planeNormal)).data());
  
  //array_t coeffs = (array_t)map1.colwise().sum()/(array_t)map2.colwise().sum();
  
  // Loopless
  //Eigen::Map<matrix_t> rayMap (rayVector.data());
  
  //rayVector.template segment<kDIM/3>(0) = std::move(((array_t)rayMap.row(0)).cwiseProduct(coeffs));
  //rayVector.template segment<kDIM/3>(kDIM/3) = std::move(((array_t)rayMap.row(1)).cwiseProduct(coeffs));
  //rayVector.template segment<kDIM/3>(2*kDIM/3) = std::move(((array_t)rayMap.row(2)).cwiseProduct(coeffs));
  
  // Semi loopless
  //for (unsigned i = 0; i < kDIM; i += 3) {
  //  rayVector.block(3,1,i,0) *= coeffs(i/3, 1);
  //}
  typedef Eigen::Array <scalar_t, kDIM, 1> array_t;
  
  array_t dat_arr1 = (rayPoint - planePoint).cwiseProduct(planeNormal);
  array_t dat_arr2 = rayVector.cwiseProduct(planeNormal);
  scalar_t sums1[kDIM], sums2[kDIM];
  for (int i = 0; i < kDIM; i += 1) {
    sums1[i] = dat_arr1[i] + dat_arr1[i+1] + dat_arr1[i+2];
    sums2[i] = dat_arr2[i] + dat_arr2[i+1] + dat_arr2[i+2];;
  }
  for (int i = 0; i < kDIM; i += 1) {
    sums1[i+1] = sums1[i+1] + dat_arr1[i] - dat_arr1[i+3];
    sums2[i+1] = sums2[i+1] + dat_arr2[i] - dat_arr2[i+3];
  }
  for (int i = 0; i < kDIM; i += 2) {
    sums1[i+2] = sums1[i+2] + dat_arr1[i] + dat_arr1[i+1] - dat_arr1[i+3] - dat_arr1[i+4];
    sums2[i+2] = sums2[i+2] + dat_arr2[i] + dat_arr2[i+1] - dat_arr1[i+3] - dat_arr2[i+4];
  }
  for (int i = 0; i < kDIM; i+=1) {
    rayVector[i] *= sums1[i]/sums2[i];
  }
  
  return rayPoint - rayVector;
}

BOOST_AUTO_TEST_SUITE(VectIntersect)

unsigned int tests = 1000000;

using Scalar = float;
using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

// Same starting position
Vector3 rv = Vector3(0.0, -1.0, -1.0);
Vector3 rp = Vector3(0.0, 0.0, 10.0);

// For the moment same normal vectors
Vector3 pn = Vector3(0.0, 0.0, 1.0);

// 8 planes
Vector3 pp0 = Vector3(0.0, 0.0, 5.0);
Vector3 pp1 = Vector3(0.0, 0.0, 6.0);
Vector3 pp2 = Vector3(0.0, 0.0, 7.0);
Vector3 pp3 = Vector3(0.0, 0.0, 8.0);
Vector3 pp4 = Vector3(0.0, 0.0, 9.0);
Vector3 pp5 = Vector3(0.0, 0.0, 10.0);
Vector3 pp6 = Vector3(0.0, 0.0, 11.0);
Vector3 pp7 = Vector3(0.0, 0.0, 12.0);
Vector3 pp8 = Vector3(0.0, 0.0, 13.0);

std::array<Vector3, 9> pps = {pp0, pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8};

template <unsigned int kPlanes> void intersectSingle() {
  for (unsigned int nt = 0; nt < tests; ++nt) {
    for (unsigned int ip = 0; ip < kPlanes; ++ip) {
      auto ips = intersect<Scalar, 3>(rv, rp, pn, pps[ip]);
    }
  }
}

template <unsigned int kPlanes> void intersectMultiple() {

  constexpr unsigned int kFullDim = kPlanes * 3;

  using VectorM = Eigen::Matrix<Scalar, kFullDim, 1>;

  VectorM rvM;
  VectorM rpM;
  VectorM pnM;
  VectorM ppM;

  for (unsigned int ip = 0; ip < kPlanes; ++ip) {
    rvM.template segment<3>(ip * 3) = rv;
    rpM.template segment<3>(ip * 3) = rp;
    pnM.template segment<3>(ip * 3) = pn;
    ppM.template segment<3>(ip * 3) = pps[ip];
  }
  for (unsigned int nt = 0; nt < tests; ++nt) {
    auto ipM = intersect<Scalar, kFullDim>(rvM, rpM, pnM, ppM);
  }
}

BOOST_AUTO_TEST_CASE(SingleIntersection4) { intersectSingle<4>(); }

BOOST_AUTO_TEST_CASE(MultiIntersection4) { intersectMultiple<4>(); }

BOOST_AUTO_TEST_CASE(SingleIntersection6) { intersectSingle<6>(); }

BOOST_AUTO_TEST_CASE(MultiIntersection6) { intersectMultiple<6>(); }

BOOST_AUTO_TEST_CASE(SingleIntersection8) { intersectSingle<8>(); }

BOOST_AUTO_TEST_CASE(MultiIntersection8) { intersectMultiple<8>(); }

BOOST_AUTO_TEST_CASE(SingleIntersection9) { intersectSingle<9>(); }

BOOST_AUTO_TEST_CASE(MultiIntersection9) { intersectMultiple<9>(); }


BOOST_AUTO_TEST_SUITE_END()

} // namespace Test

} // namespace vecint
