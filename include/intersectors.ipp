 /**
 * author: joana.niermann@cern.ch
 **/
namespace vec_intr {


auto eig_intersect_4D(ray_data<Vector4_s>& ray,
                      Vector4_s& planeNormal,
                      Vector4_s& planePoint) {
  Scalar nom   = (ray.point - planePoint).dot(planeNormal);
  Scalar coeff = nom / ray.direction.dot(planeNormal);

  intersection<Scalar, Vector4_s> results = {.path = Vector4_s(ray.point - coeff*ray.direction),
                                             .dist = coeff};       
  return std::move(results);
}



template<typename vector_t>
void eig_intersect_4D(ray_data<Vector4_s>& ray,
                      plane_data<aligned::vector<vector_t>> planes,
                      aligned::vector<intersection<Scalar, Vector4_s>> &results) {
  #ifdef DEBUG
  //TODO: make constexpr
  if (planePoints.size() != planeNormals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  for (size_t i = 0; i < planes.points.size(); i++) {
    Scalar nom   = (ray.point - planes.points[i].obj).dot(planes.normals[i].obj);
    Scalar coeff = nom / ray.direction.dot(planes.normals[i].obj);

    results.emplace_back(intersection<Scalar, Vector4_s>{.path = Vector4_s(ray.point - coeff*ray.direction), .dist = coeff});
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
void vc_intersect_vert(ray_data<Vector4_s> &ray,
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

    auto denom = ray_dir * plane_normal;
    auto nom   = (ray_point - plane_point) * plane_normal;

    scalar_t coeff = nom.sum() / denom.sum();
    auto path = (ray_point - coeff*ray_dir);

    results.emplace_back(intersection<scalar_t, simd_vec_t>{.path = std::move(path), .dist = coeff});
  }
}


template<typename vector_v, typename vector_s>
auto vc_intersect_vert(ray_data<Vector4_s> &ray,
                       Vector4_s& planeNormal,
                       Vector4_s& planePoint) {
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

  simd_vec_t  plane_normal = simd_vec_t(planeNormal.data());
  simd_vec_t  plane_point  = simd_vec_t(planePoint.data());

  auto denom = ray_dir * plane_normal;
  auto nom   = (ray_point - plane_point) * plane_normal;

  scalar_t coeff = nom.sum() / denom.sum();
  auto path = (ray_point - coeff*ray_dir);

  intersection<scalar_t, simd_vec_t> results = {.path = std::move(path), .dist = coeff};
  return std::move(results);
}



template<typename scalar_v>
void vc_intersect_hybrid(ray_data<Vector3<scalar_v>>& ray,
                         plane_data<aligned::vector<Vector3<typename scalar_v::value_type>>> &planes,
                         aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results) {

  using scalar_t = typename scalar_v::value_type;

  #ifdef DEBUG
  //TODO: make constexpr
  if (pns_struct.size() != pps_struct.size()) {
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

    scalar_v denom_x (ray.direction.x * pns_x);
    scalar_v denom_y (ray.direction.y * pns_y);
    scalar_v denom_z (ray.direction.z * pns_z);

    scalar_v denoms (denom_x + denom_y + denom_z);

    scalar_v nom_x ((ray.point.x - pps_x) * pns_x);
    scalar_v nom_y ((ray.point.y - pps_y) * pns_y);
    scalar_v nom_z ((ray.point.z - pps_z) * pns_z);

    scalar_v coeffs ((nom_x + nom_y + nom_z) / denoms);

    Vector3<scalar_v> path = {.x = (ray.point.x - coeffs*ray.direction.x),
                              .y = (ray.point.y - coeffs*ray.direction.y),
                              .z = (ray.point.z - coeffs*ray.direction.z)};

    results.emplace_back(intersection<scalar_v, Vector3<scalar_v>>{.path = std::move(path), .dist = std::move(coeffs)});
  }
}


template<typename scalar_v>
auto vc_intersect_hybrid(ray_data<Vector3<scalar_v>>& ray,
                         plane_data<Vector3<scalar_v>>& planes) {

  #ifdef DEBUG
  //TODO: make constexpr
  if (pns_struct.size() != pps_struct.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  scalar_v denom_x (ray.direction.x * planes.normals.x);
  scalar_v denom_y (ray.direction.y * planes.normals.y);
  scalar_v denom_z (ray.direction.z * planes.normals.z);

  scalar_v denoms (denom_x + denom_y + denom_z);

  scalar_v nom_x ((ray.point.x - planes.points.x) * planes.normals.x);
  scalar_v nom_y ((ray.point.y - planes.points.y) * planes.normals.y);
  scalar_v nom_z ((ray.point.z - planes.points.z) * planes.normals.z);

  scalar_v coeffs ((nom_x + nom_y + nom_z) / denoms);

  Vector3<scalar_v> path = {.x = (ray.point.x - coeffs*ray.direction.x),
                            .y = (ray.point.y - coeffs*ray.direction.y),
                            .z = (ray.point.z - coeffs*ray.direction.z)};

  intersection<scalar_v, Vector3<scalar_v>> results = {.path = std::move(path), .dist = std::move(coeffs)};
  return std::move(results);
}


template<typename scalar_v, typename vector_s, size_t kDIM = 3>
void vc_intersect_horiz(ray_data<Vector3<scalar_v>>& ray,
                        plane_data<aligned::vector<MatrixV<vector_s>>>& planes,
                        aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results,
                        size_t padding = 0) {
  #ifdef DEBUG
  //TODO: make constexpr
  if (planePoints.size() != planeNormals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  // Access to raw data that will be loaded as scalar_v
  size_t n_bytes = planes.points.size() * (planes.points.front().n_elemts() + padding);
  size_t offset  = scalar_v::Size + padding;
  if (n_bytes % (kDIM*offset) != 0) std::cout << "Warning: Input container size is not a multiple simd vector size." << std::endl;
  size_t n_inters = n_bytes / (kDIM*offset);
  if (results.capacity() < n_inters) results.reserve(n_inters);

  auto pp_ptr = planes.points.front().data();
  auto pn_ptr = planes.normals.front().data();
  scalar_v planePoint (pp_ptr);
  scalar_v planeNormal(pn_ptr);
  for (size_t i = 0; i < n_inters; i++) {

    scalar_v denoms(ray.direction.x * planeNormal);
    scalar_v coeffs((ray.point.x - planePoint) * planeNormal);
    planeNormal.load(pn_ptr += offset, Vc::Streaming);
    planePoint.load( pp_ptr += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.y, planeNormal, denoms);
    coeffs = Vc::fma((ray.point.y - planePoint), planeNormal, coeffs);
    planeNormal.load(pn_ptr += offset, Vc::Streaming);
    planePoint.load( pp_ptr += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.z, planeNormal, denoms);
    coeffs = Vc::fma((ray.point.z - planePoint), planeNormal, coeffs);
    planeNormal.load(pn_ptr += offset, Vc::Streaming);
    planePoint.load( pp_ptr += offset, Vc::Streaming);
    coeffs /= denoms;

    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, -ray.direction.x, ray.point.x), 
                              .y = Vc::fma(coeffs, -ray.direction.y, ray.point.y), 
                              .z = Vc::fma(coeffs, -ray.direction.z, ray.point.z)};

    results.emplace_back(intersection<scalar_v, Vector3<scalar_v>>{.path = std::move(path), .dist = std::move(coeffs)});
  }
}


template<typename scalar_v, typename data_ptr_t, size_t kDIM = 3>
auto vc_intersect_horiz(ray_data<Vector3<scalar_v>>& ray,
                        data_ptr_t pl_normals_ptr,
                        data_ptr_t pl_points_ptr,
                        size_t& offset) {
    #ifdef DEBUG
    if (pl_normals_ptr == nullptr || pl_points_ptr == nullptr) {
      std::cerr << "Passed invalid data collection pointer to intersection" << std::endl;
      return intersection<scalar_v, Vector3<scalar_v>>{};
    }
    #endif

    scalar_v planePoint (pl_points_ptr);
    scalar_v planeNormal(pl_normals_ptr);

    scalar_v denoms(ray.direction.x * planeNormal);
    scalar_v coeffs((ray.point.x - planePoint) * planeNormal);
    planeNormal.load(pl_normals_ptr += offset, Vc::Streaming);
    planePoint.load( pl_points_ptr  += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.y, planeNormal, denoms);
    coeffs = Vc::fma((ray.point.y - planePoint), planeNormal, coeffs);
    planeNormal.load(pl_normals_ptr += offset, Vc::Streaming);
    planePoint.load( pl_points_ptr  += offset, Vc::Streaming);

    denoms = Vc::fma(ray.direction.z, planeNormal, denoms);
    coeffs = Vc::fma((ray.point.z - planePoint), planeNormal, coeffs);
    planeNormal.load(pl_normals_ptr += offset, Vc::Streaming);
    planePoint.load( pl_points_ptr  += offset, Vc::Streaming);
    coeffs /= denoms;

    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, -ray.direction.x, ray.point.x), 
                              .y = Vc::fma(coeffs, -ray.direction.y, ray.point.y), 
                              .z = Vc::fma(coeffs, -ray.direction.z, ray.point.z)};

    intersection<scalar_v, Vector3<scalar_v>> results = {.path = std::move(path), .dist = std::move(coeffs)};
    return std::move(results);
}


} // namespace vec_intr
