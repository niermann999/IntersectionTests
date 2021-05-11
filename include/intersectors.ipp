 /**
 * author: joana.niermann@cern.ch
 **/
namespace vec_intr {


template<typename vector_s>
inline auto eig_intersect_4D(ray_data<vector_s>& ray,
                             plane_data<vector_s> &plane) {

  using scalar_t = typename vector_s::scalar_type;
  using vector_t = typename vector_s::type;

  Scalar coeff = (plane.points() - ray.point()).dot(plane.normals());
  Scalar denom = ray.direction().dot(plane.normals());
  if (denom == 0) return intersection<scalar_t, vector_t>{};
  coeff /= denom;

  intersection<scalar_t, vector_t> results = {.path = vector_t(ray.point() + coeff*ray.direction()),
                                              .dist = coeff};       
  return std::move(results);
}


template<typename vector_s>
inline void eig_intersect_4D(ray_data<vector_s> &ray,
                             aligned::vector<plane_data<vector_s> > &planes,
                             aligned::vector<intersection<typename vector_s::scalar_type, typename vector_s::type> > &results) {

  using scalar_t = typename vector_s::scalar_type;
  using vector_t = typename vector_s::type;

  for (auto &plane : planes) {
    scalar_t coeff = (plane.points()- ray.point()).dot(plane.normals());
    scalar_t denom = ray.direction().dot(plane.normals());
    if (denom == 0) return;
    coeff /= denom;

    results.emplace_back(intersection<scalar_t, vector_t>{.path = vector_t(ray.point() + coeff*ray.direction()), .dist = coeff});
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

template<typename vector_s>
inline void vc_intersect_vert(ray_data<vector_s> &ray,
                              aligned::vector<plane_data<vector_s> > &planes,
                              aligned::vector<intersection<typename vector_s::scalar_type, Vc::SimdArray<typename vector_s::scalar_type, 4> > > &results) {

  using scalar_t = typename vector_s::scalar_type;
  using simd_vec_t = Vc::SimdArray<scalar_t, 4>;

  auto ray_dir   = simd_vec_t(ray.direction.data());
  auto ray_point = simd_vec_t(ray.point.data());

  simd_vec_t plane_normal, plane_point;

  for (auto &plane : planes) {
    plane_normal = simd_vec_t(plane.normals.data());
    plane_point  = simd_vec_t(plane.points.data());

    scalar_t coeff = ((plane_point - ray_point) * plane_normal).sum();
    scalar_t denom = (ray_dir * plane_normal).sum();
    if (denom == 0) return;
    coeff /= denom;

    auto path = (ray_point + coeff*ray_dir);

    results.emplace_back(intersection<scalar_t, simd_vec_t>{.path = path, .dist = coeff});
  }
}


template<typename vector_s>
inline auto vc_intersect_vert(ray_data<vector_s> &ray,
                              plane_data<vector_s> &plane) {

  using scalar_t = typename vector_s::scalar_type;
  using simd_vec_t = Vc::SimdArray<scalar_t, 4>;

  auto ray_dir   = simd_vec_t(ray.direction.data());
  auto ray_point = simd_vec_t(ray.point.data());

  simd_vec_t  plane_normal = simd_vec_t(plane.normals.data());
  simd_vec_t  plane_point  = simd_vec_t(plane.points.data());

  scalar_t coeff = ((plane_point - ray_point) * plane_normal).sum();
  scalar_t denom = (ray_dir * plane_normal).sum();
  if (denom == 0) return intersection<scalar_t, simd_vec_t>{};
  coeff /= denom;

  auto path = (ray_point + coeff*ray_dir);

  intersection<scalar_t, simd_vec_t> results = {.path = path, .dist = coeff};
  return std::move(results);
}


/*template<typename vector_v>
inline void vc_intersect_hybrid(ray_data<vector_v> &ray,
                                plane_data<aligned::vector<Vector3<typename vector_v::scalar_type> > > &planes,
                                std::array<intersection<typename vector_v::vec_type, typename vector_v::type>, 10> &results) {*/
template<typename vector_v>
inline void vc_intersect_hybrid(ray_data<vector_v> &ray,
                                plane_data<aligned::vector<Vector3<typename vector_v::scalar_type> > > &planes,
                                aligned::vector<intersection<typename vector_v::vec_type, typename vector_v::type> > &results) {

  using scalar_t = typename vector_v::scalar_type;
  using scalar_v = typename vector_v::vec_type;
  using output_t = typename vector_v::type;

  #ifdef DEBUG
  //TODO: make constexpr
  if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }
  #endif

  //size_t i = 0;
  // Vector iterate
  for (Index_v i(scalar_v::IndexesFromZero()); (i < Index_v(planes.points.size())).isFull(); i += Index_v(scalar_v::Size)) {
    scalar_v pps_x = planes.points[i][&Vector3<scalar_t>::x];
    scalar_v pps_y = planes.points[i][&Vector3<scalar_t>::y];
    scalar_v pps_z = planes.points[i][&Vector3<scalar_t>::z];

    scalar_v pns_x = planes.normals[i][&Vector3<scalar_t>::x];
    scalar_v pns_y = planes.normals[i][&Vector3<scalar_t>::y];
    scalar_v pns_z = planes.normals[i][&Vector3<scalar_t>::z];

    scalar_v denoms (ray.direction().x * pns_x);
    scalar_v coeffs ((pps_x - ray.point().x) * pns_x);

    denoms = Vc::fma(ray.direction().y, pns_y, denoms);
    coeffs = Vc::fma((pps_y - ray.point().y), pns_y, coeffs);

    denoms = Vc::fma(ray.direction().z, pns_z, denoms);
    coeffs = Vc::fma((pps_z - ray.point().z), pns_z, coeffs);
    coeffs /= denoms;

    auto check_sum = coeffs.sum();
    if (std::isnan(check_sum) || std::isinf(check_sum)) return;

    output_t path = {.x = Vc::fma(coeffs, ray.direction().x, ray.point().x), 
                                    .y = Vc::fma(coeffs, ray.direction().y, ray.point().y), 
                                    .z = Vc::fma(coeffs, ray.direction().z, ray.point().z)};

    //results.emplace_back(intersection<scalar_v, output_t>{.path = path, .dist = coeffs});
    intersection<scalar_v, output_t> bla = {.path = path, .dist = coeffs};
    results[0] = bla;
    //i += 1;
  }
}


template<typename vector_v>
inline auto vc_intersect_hybrid(ray_data<vector_v> &ray,
                                plane_data<vector_v> &planes) {
                                    
  using scalar_v = typename vector_v::vec_type;
  using output_t = typename vector_v::type;

  scalar_v denoms (ray.direction().x * planes.normals().x);
  scalar_v coeffs ((planes.points().x - ray.point().x) * planes.normals().x);

  denoms = Vc::fma(ray.direction().y, planes.normals().y, denoms);
  coeffs = Vc::fma((planes.points().y - ray.point().y), planes.normals().y, coeffs);

  denoms = Vc::fma(ray.direction().z, planes.normals().z, denoms);
  coeffs = Vc::fma((planes.points().z - ray.point().z), planes.normals().z, coeffs);
  coeffs /= denoms;

  auto check_sum = coeffs.sum();
  if (std::isnan(check_sum) || std::isinf(check_sum)) return intersection<scalar_v, Vector3<scalar_v> >{};

  Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction().x, ray.point().x), 
                            .y = Vc::fma(coeffs, ray.direction().y, ray.point().y), 
                            .z = Vc::fma(coeffs, ray.direction().z, ray.point().z)};

  intersection<scalar_v, Vector3<scalar_v> > results = {.path = path, .dist = coeffs};
  return std::move(results);
}


template<typename vector_v>
inline void vc_intersect_horiz(ray_data<vector_v> &ray,
                               aligned::vector<plane_data<vector_v> > &planes,
                               aligned::vector<intersection<typename vector_v::vec_type, typename vector_v::type> > &results) {

  using scalar_v = typename vector_v::vec_type;
  using output_t = typename vector_v::type;

  for (auto &plane : planes) {
    scalar_v denoms(ray.direction().x * plane.normals().x);
    scalar_v coeffs((plane.points().x - ray.point().x) * plane.normals().x);

    denoms = Vc::fma(ray.direction().y, plane.normals().y, denoms);
    coeffs = Vc::fma((plane.points().y - ray.point().y), plane.normals().y, coeffs);

    denoms = Vc::fma(ray.direction().z, plane.normals().z, denoms);
    coeffs = Vc::fma((plane.points().z - ray.point().z), plane.normals().z, coeffs);
    coeffs /= denoms;

    auto check_sum = coeffs.sum();
    if (std::isnan(check_sum) || std::isinf(check_sum)) return;


    output_t path = {.x = Vc::fma(coeffs, ray.direction().x, ray.point().x), 
                     .y = Vc::fma(coeffs, ray.direction().y, ray.point().y), 
                     .z = Vc::fma(coeffs, ray.direction().z, ray.point().z)};

    results.emplace_back(intersection<scalar_v, output_t>{.path = path, .dist = coeffs});
  }
}


template<typename vector_v>
inline auto vc_intersect_horiz(ray_data<vector_v> &ray,
                               plane_data<vector_v> &plane) {

    using scalar_v = typename vector_v::vec_type;
    using output_t = typename vector_v::type;

    scalar_v denoms(ray.direction().x * plane.normals().x);
    scalar_v coeffs((plane.points().x - ray.point().x) * plane.normals().x);

    denoms = Vc::fma(ray.direction().y, plane.normals().y, denoms);
    coeffs = Vc::fma((plane.points().y - ray.point().y), plane.normals().y, coeffs);

    denoms = Vc::fma(ray.direction().z, plane.normals().z, denoms);
    coeffs = Vc::fma((plane.points().z - ray.point().z), plane.normals().z, coeffs);
    coeffs /= denoms;
    
    auto check_sum = coeffs.sum();
    if (std::isnan(check_sum) || std::isinf(check_sum)) {
      return intersection<Scalar_v, typename vector_v::type> {};
    }

    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction().x, ray.point().x), 
                              .y = Vc::fma(coeffs, ray.direction().y, ray.point().y), 
                              .z = Vc::fma(coeffs, ray.direction().z, ray.point().z)};

    intersection<scalar_v, output_t> results = {.path = path, .dist = coeffs};
    return std::move(results);
}

} // namespace vec_intr
