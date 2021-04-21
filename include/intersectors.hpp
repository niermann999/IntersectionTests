#include <types.hpp>

#pragma once

namespace vec_intr {
//-----------------------------------------------------------------------
// Eigen
//-----------------------------------------------------------------------
// Naive Eigen implementation
auto eig_intersect_4D(ray_data<Vector4_s>& ray,
                      Vector4_s& planeNormal,
                      Vector4_s& planePoint) ;


// Process all intersections at once
template<typename vector_t>
void eig_intersect_4D(ray_data<Vector4_s>& ray,
                      plane_data<aligned::vector<vector_t>> planes,
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
void vc_intersect_vert(ray_data<Vector4_s> &ray,
                       plane_data<aligned::vector<vector_s>> &planes,
                       aligned::vector<intersection<typename vector_v::value_type, Vc::SimdArray<typename vector_v::value_type, 4>>> &results);


template<typename vector_v, typename vector_s>
auto vc_intersect_vert(ray_data<Vector4_s> &ray,
                       Vector4_s& planeNormal,
                       Vector4_s& planePoint);


// Vertical vectorization on data as it is, using interleaved memory wrapper
template<typename scalar_v>
void vc_intersect_hybrid(ray_data<Vector3<scalar_v>>& ray,
                         plane_data<aligned::vector<Vector3<typename scalar_v::value_type>>> &planes,
                         aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results);


template<typename scalar_v>
auto vc_intersect_hybrid(ray_data<Vector3<scalar_v>>& ray,
                         plane_data<Vector3<scalar_v>>& planes);


// This is terrible!
/*template<typename vector_v, typename vector_s> // Set type traits according to assumptions that are made
void vc_intersect_vert(//ray_data<vector_s> &ray,
                       ray_data<Vector4_s> &ray,
                       //mem_t<kDIM> &planePoints) {    // needs explicit copy, because fromRawData() does not work (uses 'offsetof' on private member)
                       //Vc::span<vector_t, kDIM> &planePoints) { // terribly slow somehow
                       plane_data<aligned::vector<vector_s>> &planes,
                       aligned::vector<intersection<std::vector<typename vector_v::value_type>, aligned::vector<vector_v>>> &results,
                       size_t stride = vector_v::Size);
*/


template<typename scalar_v>
void vc_intersect_hybrid(Vector3<scalar_v> &rayVector,
                         Vector3<scalar_v> &rayPoint,
                         aligned::vector<Vector3<typename scalar_v::value_type>> &pns_struct,
                         aligned::vector<Vector3<typename scalar_v::value_type>> &pps_struct,
                         aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results); 


template<typename scalar_v, typename vector_s, size_t kDIM = 4>
void vc_intersect_horiz(ray_data<Vector3<scalar_v>>& ray,
                        plane_data<aligned::vector<MatrixV<vector_s>>>& planes,
                        aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results); 


template<typename scalar_v, typename data_ptr_t, size_t kDIM = 3>
auto vc_intersect_horiz(ray_data<Vector3<scalar_v>>& ray,
                        data_ptr_t pp_ptr,
                        data_ptr_t pn_ptr,
                        size_t& offset); 


template<typename scalar_v,typename vector_v, size_t kDIM = 3>
bool initialize_data_ptrs (Scalar **plane_normals,
                           Scalar **plane_points,
                           aligned::vector<vector_v>& normals_vec,
                           aligned::vector<vector_v>& points_vec); 
} //namespace vec_intr

#include <intersectors.ipp>