 /**
 * author: joana.niermann@cern.ch
 **/
namespace vec_intr {


inline auto eig_intersect_4D(ray_data<Vector4_s> &ray,
                             Vector4_s &planeNormal,
                             Vector4_s &planePoint) {
  Scalar coeff = (planePoint - ray.point).dot(planeNormal);
  Scalar denom = ray.direction.dot(planeNormal);
  if (denom == 0) return intersection<Scalar, Vector4_s>{};
  coeff /= denom;

  intersection<Scalar, Vector4_s> results = {.path = Vector4_s(ray.point + coeff*ray.direction),
                                             .dist = coeff};       
  return std::move(results);
}


template<typename vector_t>
inline void eig_intersect_4D(ray_data<Vector4_s> &ray,
                             plane_data<aligned::vector<vector_t>> &planes,
                             aligned::vector<intersection<Scalar, Vector4_s>> &results) {
  #ifdef DEBUG
  //TODO: make constexpr
  if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  for (size_t i = 0; i < planes.points.size(); i++) {
    Scalar coeff = (planes.points[i].obj - ray.point).dot(planes.normals[i].obj);
    Scalar denom = ray.direction.dot(planes.normals[i].obj);
    if (denom == 0) return;
    coeff /= denom;

    results.emplace_back(intersection<Scalar, Vector4_s>{.path = Vector4_s(ray.point + coeff*ray.direction), .dist = coeff});
  }
}



/*template <typename scalar_t, size_t kDIM>
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
}*/


//-----------------------------------------------------------------------
// Vc
//-----------------------------------------------------------------------

template<typename vector_v, typename vector_s>
inline void vc_intersect_vert(ray_data<Vector4_s> &ray,
                              plane_data<aligned::vector<vector_s>> &planes,
                              aligned::vector<intersection<typename vector_v::value_type, Vc::SimdArray<typename vector_v::value_type, 4>>> &results) {
  using scalar_t = typename vector_v::value_type;
  using simd_vec_t = Vc::SimdArray<scalar_t, 4>;

  #ifdef DEBUG
  //TODO: make constexpr
  if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  auto ray_dir   = simd_vec_t(ray.direction.data());
  auto ray_point = simd_vec_t(ray.point.data());

  simd_vec_t plane_normal, plane_point;

  for (size_t i = 0; i < planes.normals.size(); i++) {
    plane_normal = simd_vec_t(planes.normals[i].data());
    plane_point  = simd_vec_t(planes.points[i].data());

    scalar_t coeff = ((plane_point - ray_point) * plane_normal).sum();
    scalar_t denom = (ray_dir * plane_normal).sum();
    if (denom == 0) return;
    coeff /= denom;

    auto path = (ray_point + coeff*ray_dir);

    results.emplace_back(intersection<scalar_t, simd_vec_t>{.path = path, .dist = coeff});
  }
}


template<typename vector_v, typename vector_s>
inline auto vc_intersect_vert(ray_data<Vector4_s> &ray,
                              Vector4_s &planeNormal,
                              Vector4_s &planePoint) {
  using scalar_t = typename vector_v::value_type;
  using simd_vec_t = Vc::SimdArray<scalar_t, 4>;

  auto ray_dir   = simd_vec_t(ray.direction.data());
  auto ray_point = simd_vec_t(ray.point.data());

  simd_vec_t  plane_normal = simd_vec_t(planeNormal.data());
  simd_vec_t  plane_point  = simd_vec_t(planePoint.data());

  scalar_t coeff = ((plane_point - ray_point) * plane_normal).sum();
  scalar_t denom = (ray_dir * plane_normal).sum();
  if (denom == 0) return intersection<scalar_t, simd_vec_t>{};
  coeff /= denom;

  auto path = (ray_point + coeff*ray_dir);

  intersection<scalar_t, simd_vec_t> results = {.path = path, .dist = coeff};
  return std::move(results);
}



template<typename scalar_v>
inline void vc_intersect_hybrid(ray_data<Vector3<scalar_v>> &ray,
                                plane_data<aligned::vector<Vector3<typename scalar_v::value_type>>> &planes,
                                aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results) {

  using scalar_t = typename scalar_v::value_type;

  #ifdef DEBUG
  //TODO: make constexpr
  if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  // Vector iterate
  for (Index_v i(scalar_v::IndexesFromZero()); (i < Index_v(planes.points.size())).isFull(); i += Index_v(scalar_v::Size)) {
    scalar_v pps_x = planes.points[i][&Vector3<scalar_t>::x];
    scalar_v pps_y = planes.points[i][&Vector3<scalar_t>::y];
    scalar_v pps_z = planes.points[i][&Vector3<scalar_t>::z];

    scalar_v pns_x = planes.normals[i][&Vector3<scalar_t>::x];
    scalar_v pns_y = planes.normals[i][&Vector3<scalar_t>::y];
    scalar_v pns_z = planes.normals[i][&Vector3<scalar_t>::z];

    scalar_v denoms (ray.direction.x * pns_x);
    scalar_v coeffs ((pps_x - ray.point.x) * pns_x);

    denoms = Vc::fma(ray.direction.y, pns_y, denoms);
    coeffs = Vc::fma((pps_y - ray.point.y), pns_y, coeffs);

    denoms = Vc::fma(ray.direction.z, pns_z, denoms);
    coeffs = Vc::fma((pps_z - ray.point.z), pns_z, coeffs);
    coeffs /= denoms;

    auto check_sum = coeffs.sum();
    if (std::isnan(check_sum) || std::isinf(check_sum)) return;

    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction.x, ray.point.x), 
                              .y = Vc::fma(coeffs, ray.direction.y, ray.point.y), 
                              .z = Vc::fma(coeffs, ray.direction.z, ray.point.z)};

    results.emplace_back(intersection<scalar_v, Vector3<scalar_v>>{.path = path, .dist = coeffs});
  }
}


template<typename scalar_v>
inline auto vc_intersect_hybrid(ray_data<Vector3<scalar_v>> &ray,
                                plane_data<Vector3<scalar_v>> &planes) {

  scalar_v denoms (ray.direction.x * planes.normals.x);
  scalar_v coeffs ((planes.points.x - ray.point.x) * planes.normals.x);

  denoms = Vc::fma(ray.direction.y, planes.normals.y, denoms);
  coeffs = Vc::fma((planes.points.y - ray.point.y), planes.normals.y, coeffs);

  denoms = Vc::fma(ray.direction.z, planes.normals.z, denoms);
  coeffs = Vc::fma((planes.points.z - ray.point.z), planes.normals.z, coeffs);
  coeffs /= denoms;

  auto check_sum = coeffs.sum();
  if (std::isnan(check_sum) || std::isinf(check_sum)) return intersection<scalar_v, Vector3<scalar_v>>{};

  Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction.x, ray.point.x), 
                            .y = Vc::fma(coeffs, ray.direction.y, ray.point.y), 
                            .z = Vc::fma(coeffs, ray.direction.z, ray.point.z)};

  intersection<scalar_v, Vector3<scalar_v>> results = {.path = path, .dist = coeffs};
  return std::move(results);
}

template<typename scalar_v, typename vector_s, size_t kDIM = 3>
inline void vc_intersect_horiz(ray_data<Vector3<scalar_v>> &ray,
                               plane_data<aligned::vector<MatrixV<vector_s>>> &planes,
                               aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results,
                               size_t padding = 0) {
  #ifdef DEBUG
  //TODO: make constexpr
  if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  // Access to raw data that will be loaded as scalar_v
  size_t n_float_pnt = planes.points.size() * (planes.points.front().n_elemts() + padding);
  size_t offset  = planes.points.front().n_elemts() + padding;
  // Process 3 geometrical coordinates
  if (n_float_pnt % (kDIM*offset) != 0) std::cout << "Warning: Input container size is not a multiple simd vector size." << std::endl;
  size_t n_loops = n_float_pnt / (kDIM*offset);
  if (results.capacity() < n_loops) results.reserve(n_loops);

  auto pp_ptr = planes.points.front().data();
  auto pn_ptr = planes.normals.front().data();
  scalar_v planePoint (pp_ptr);
  scalar_v planeNormal(pn_ptr);
  for (size_t i = 0; i < n_loops; i++) {

    scalar_v denoms(ray.direction.x * planeNormal);
    scalar_v coeffs((planePoint - ray.point.x) * planeNormal);
    planeNormal.load(pn_ptr += offset, Vc::Streaming);
    planePoint.load( pp_ptr += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.y, planeNormal, denoms);
    coeffs = Vc::fma((planePoint - ray.point.y), planeNormal, coeffs);
    planeNormal.load(pn_ptr += offset, Vc::Streaming);
    planePoint.load( pp_ptr += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.z, planeNormal, denoms);
    coeffs = Vc::fma((planePoint - ray.point.z), planeNormal, coeffs);
    planeNormal.load(pn_ptr += offset, Vc::Streaming);
    planePoint.load( pp_ptr += offset, Vc::Streaming);
    coeffs /= denoms;

    auto check_sum = coeffs.sum();
    if (std::isnan(check_sum) || std::isinf(check_sum)) return;


    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction.x, ray.point.x), 
                              .y = Vc::fma(coeffs, ray.direction.y, ray.point.y), 
                              .z = Vc::fma(coeffs, ray.direction.z, ray.point.z)};

    results.emplace_back(intersection<scalar_v, Vector3<scalar_v>>{.path = path, .dist = coeffs});
  }
}


template<typename scalar_v, typename data_ptr_t, size_t kDIM = 3>
inline auto vc_intersect_horiz(ray_data<Vector3<scalar_v>> &ray,
                               data_ptr_t pl_normals_ptr,
                               data_ptr_t pl_points_ptr,
                               size_t offset) {
    #ifdef DEBUG
    if (pl_normals_ptr == nullptr || pl_points_ptr == nullptr) {
      std::cerr << "Passed invalid data collection pointer to intersection" << std::endl;
      return intersection<scalar_v, Vector3<scalar_v>>{};
    }
    #endif

    scalar_v planePoint (pl_points_ptr);
    scalar_v planeNormal(pl_normals_ptr);

    scalar_v denoms(ray.direction.x * planeNormal);
    scalar_v coeffs((planePoint - ray.point.x) * planeNormal);
    planeNormal.load(pl_normals_ptr += offset, Vc::Streaming);
    planePoint.load( pl_points_ptr  += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.y, planeNormal, denoms);
    coeffs = Vc::fma((planePoint - ray.point.y), planeNormal, coeffs);
    planeNormal.load(pl_normals_ptr += offset, Vc::Streaming);
    planePoint.load( pl_points_ptr  += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.z, planeNormal, denoms);
    coeffs = Vc::fma((planePoint - ray.point.z), planeNormal, coeffs);
    planeNormal.load(pl_normals_ptr += offset, Vc::Streaming);
    planePoint.load( pl_points_ptr  += offset, Vc::Streaming);
    coeffs /= denoms;

    auto check_sum = coeffs.sum();
    if (std::isnan(check_sum) || std::isinf(check_sum)) {
      return intersection<scalar_v, Vector3<scalar_v>> {};
    }

    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction.x, ray.point.x), 
                              .y = Vc::fma(coeffs, ray.direction.y, ray.point.y), 
                              .z = Vc::fma(coeffs, ray.direction.z, ray.point.z)};

    intersection<scalar_v, Vector3<scalar_v>> results = {.path = path, .dist = coeffs};
    return std::move(results);
}


} // namespace vec_intr
