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
inline auto eig_intersect_4D(ray_data<Vector4_s>& ray,
                             Vector4_s& planeNormal,
                             Vector4_s& planePoint);


// Process all intersections at once
template<typename vector_t>
inline void eig_intersect_4D(ray_data<Vector4_s>& ray,
                             plane_data<aligned::vector<vector_t>> &planes,
                             aligned::vector<intersection<Scalar, Vector4_s>> &results);


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
template<typename vector_v, typename vector_s>
inline void vc_intersect_vert(ray_data<Vector4_s> &ray,
                              plane_data<aligned::vector<vector_s>> &planes,
                              aligned::vector<intersection<typename vector_v::value_type, Vc::SimdArray<typename vector_v::value_type, 4>>> &results);


template<typename vector_v, typename vector_s>
inline auto vc_intersect_vert(ray_data<Vector4_s> &ray,
                              Vector4_s &planeNormal,
                              Vector4_s &planePoint);


// Vertical vectorization on data as it is, using interleaved memory wrapper
template<typename scalar_v>
inline void vc_intersect_hybrid(ray_data<Vector3<scalar_v>> &ray,
                                plane_data<aligned::vector<Vector3<typename scalar_v::value_type>>> &planes,
                                aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results);


template<typename scalar_v>
inline auto vc_intersect_hybrid(ray_data<Vector3<scalar_v>> &ray,
                                plane_data<Vector3<scalar_v>> &planes);



template<typename scalar_v, typename vector_s, size_t kDIM = 3>
inline void vc_intersect_horiz(ray_data<Vector3<scalar_v>> &ray,
                               plane_data<aligned::vector<MatrixV<vector_s>>> &planes,
                               aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results); 


template<typename scalar_v, typename data_ptr_t, size_t kDIM = 3>
inline auto vc_intersect_horiz(ray_data<Vector3<scalar_v>> &ray,
                               data_ptr_t pl_normals_ptr,
                               data_ptr_t pl_points_ptr,
                               size_t offset); 


template<typename scalar_v,typename vector_v, size_t kDIM = 3>
inline bool initialize_data_ptrs (Scalar **plane_normals,
                                  Scalar **plane_points,
                                  aligned::vector<vector_v> &normals_vec,
                                  aligned::vector<vector_v> &points_vec); 
} //namespace vec_intr

#include <intersectors.ipp>