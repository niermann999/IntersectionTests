 /**
 * author: joana.niermann@cern.ch
 **/
#include <types.hpp>

#pragma once

namespace vec_intr {
//-----------------------------------------------------------------------
// Eigen
//-----------------------------------------------------------------------
// Naive Eigen implementation
template<typename vector_s>
inline auto eig_intersect_4D(ray_data<vector_s>& ray,
                             plane_data<vector_s> &plane);


// Process all intersections at once
template<typename vector_s>
inline void eig_intersect_4D(ray_data<vector_s> &ray,
                             aligned::vector<plane_data<vector_s> > &planes,
                             aligned::vector<intersection<typename vector_s::scalar_type, typename vector_s::type> > &results);


// Pass data in one big matrix and try to autovectorize
/*template <typename scalar_t, size_t kDIM>
auto intersect(Eigen::Matrix<scalar_t, kDIM, 3> rayVector,
               Eigen::Matrix<scalar_t, kDIM, 3> rayPoint,
               Eigen::Matrix<scalar_t, kDIM, 3> planeNormal,
               Eigen::Matrix<scalar_t, kDIM, 3> planePoint);*/

//-----------------------------------------------------------------------
// Vc
//-----------------------------------------------------------------------

// Vertical vectorization using simd array class
template<typename vector_s>
inline void vc_intersect_vert(ray_data<vector_s> &ray,
                              aligned::vector<plane_data<vector_s> > &planes,
                              aligned::vector<intersection<typename vector_s::scalar_type, Vc::SimdArray<typename vector_s::scalar_type, 4> > > &results);


template<typename vector_s>
inline auto vc_intersect_vert(ray_data<vector_s> &ray,
                              plane_data<vector_s> &planes);


// Vertical vectorization on data as it is, using interleaved memory wrapper
template<typename vector_v>
inline void vc_intersect_hybrid(ray_data<vector_v> &ray,
                                plane_data<aligned::vector<Vector3<typename vector_v::scalar_type> > > &planes,
                                aligned::vector<intersection<typename vector_v::vec_type, typename vector_v::type> > &results);


template<typename vector_v>
inline auto vc_intersect_hybrid(ray_data<vector_v> &ray,
                                plane_data<vector_v> &planes);


// Horizontal vectorization on interleaved vectors
template<typename vector_v>
inline void vc_intersect_horiz(ray_data<vector_v> &ray,
                               aligned::vector<plane_data<vector_v> > &planes,
                               aligned::vector<intersection<typename vector_v::vec_type, typename vector_v::type> > &results);


template<typename vector_v>
inline auto vc_intersect_horiz(ray_data<vector_v> &ray,
                               plane_data<vector_v> &plane);

} //namespace vec_intr

#include <intersectors.ipp>