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

template <typename scalar_t, unsigned int kDIM>
auto intersect(Eigen::Matrix<scalar_t, kDIM, 3> rayVector,
               Eigen::Matrix<scalar_t, kDIM, 3> rayPoint,
               Eigen::Matrix<scalar_t, kDIM, 3> planeNormal,
               Eigen::Matrix<scalar_t, kDIM, 3> planePoint) {
  Eigen::Matrix<scalar_t, kDIM, 3> tmpM_1 (std::move((rayPoint - planePoint).cwiseProduct(planeNormal)));
  Eigen::Matrix<scalar_t, kDIM, 3> tmpM_2 (std::move(rayVector.cwiseProduct(planeNormal)));
  Eigen::Array<scalar_t, kDIM, 1> coeffs ((tmpM_1.col(0) + tmpM_1.col(1) + tmpM_1.col(3)).array() / (tmpM_2.col(0) + tmpM_2.col(1) + tmpM_2.col(3)).array());
  
  // Broadcast coefficients onto ray-vector matrix
   rayVector.col(0).array() *= coeffs;
   rayVector.col(1).array() *= coeffs;
   rayVector.col(2).array() *= coeffs;
  
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
      auto ips = intersect<Scalar, 1>(rv, rp, pn, pps[ip]);
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
    ppM.row(ip) << pps[ip];
  }
  for (unsigned int nt = 0; nt < tests; ++nt) {
    auto ipM = intersect<Scalar, kPlanes>(rvM, rpM, pnM, ppM);
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
