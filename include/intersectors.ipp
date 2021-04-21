namespace vec_intr {


auto eig_intersect_4D(ray_data<Vector4_s>& ray,
                      Vector4_s& planeNormal,
                      Vector4_s& planePoint) {
  Scalar nom   = (ray.point - planePoint).dot(planeNormal);
  Scalar coeff = nom / ray.direction.dot(planeNormal);

  intersection<Scalar, Vector4_s> results = {.path = Vector4_s(ray.point - coeff*ray.direction),
                                             .dist = coeff};       
  return results;
}



template<typename vector_t>
void eig_intersect_4D(ray_data<Vector4_s>& ray,
                      plane_data<aligned::vector<vector_t>> planes,
                      aligned::vector<intersection<Scalar, Vector4_s>> &results) {
  /*if (planePoints.size() != planeNormals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }*/

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

  /*if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }*/

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

  /*if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }*/

  auto ray_dir   = simd_vec_t(ray.direction.data());
  auto ray_point = simd_vec_t(ray.point.data());

  simd_vec_t  plane_normal = simd_vec_t(planeNormal.data());
  simd_vec_t  plane_point  = simd_vec_t(planePoint.data());

  auto denom = ray_dir * plane_normal;
  auto nom   = (ray_point - plane_point) * plane_normal;

  scalar_t coeff = nom.sum() / denom.sum();
  auto path = (ray_point - coeff*ray_dir);

  intersection<scalar_t, simd_vec_t> res = {.path = std::move(path), .dist = coeff};
  return res;
}


// TODO: find safe(r) memory wrapper
// This is terrible!
/*template<typename vector_v, typename vector_s> // Set type traits according to assumptions that are made
void vc_intersect_vert(//ray_data<vector_s> &ray,
                       ray_data<Vector4_s> &ray,
                       //mem_t<kDIM> &planePoints) {    // needs explicit copy, because fromRawData() does not work (uses 'offsetof' on private member)
                       //Vc::span<vector_t, kDIM> &planePoints) { // terribly slow somehow
                       plane_data<aligned::vector<vector_s>> &planes,
                       aligned::vector<intersection<std::vector<typename vector_v::value_type>, aligned::vector<vector_v>>> &results,
                       size_t stride = vector_v::Size) {
  using scalar_t = typename vector_v::value_type;
  using result_t = aligned::vector<intersection<std::vector<scalar_t>, vector_v>>;
  using mask_t   = typename vector_v::MaskType;

  if (planes.points.size() != planes.normals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }

  // Access to raw data that will be reinterpreted as vector_t
  size_t n_bytes = planes.points.size() * planes.points.front().rows() * planes.points.front().cols();
  if (n_bytes % stride != 0) std::cout << "Warning: Input container size is not a multiple simd vector size." << std::endl;
  size_t n_inters = n_bytes / stride;
  if (results.capacity() < n_inters) results.reserve(n_inters);

  mask_t mask_ini = mask_t::One();
  size_t n_vec = 1;
  // Multiple geometric vectors fit into one compute vector
  if (4./vector_v::Size < 1) {
    mask_ini = vector_v::IndexesFromZero() < 4;
  }
  // In this case the number of compute vectors that holds one geometric vector is an integer number
  else {
    n_vec = 4/vector_v::Size;
  }

  size_t n_sub_loop = vector_v::Size/4 < 1 ? 0 : vector_v::Size/4 - 1;
  std::vector<scalar_t> coeff;
  aligned::vector<vector_v> ray_dirs, ray_points, result_rds;
  coeff.reserve(n_sub_loop + 1);
  ray_dirs.reserve(n_vec);
  ray_points.reserve(n_vec);
  result_rds.reserve(n_vec);
  
  auto pp_ptr = planes.points.front().data();
  auto pn_ptr = planes.normals.front().data();
  auto rd_ptr = ray.direction.data();
  auto rp_ptr = ray.point.data();
  vector_v plane_point(pp_ptr, Vc::Streaming);
  vector_v plane_normal(pn_ptr,Vc::Streaming);
  vector_v ray_direction(rd_ptr);
  vector_v ray_point(rp_ptr);
  scalar_t nom_sum = 0, denom_sum = 0;
  for (size_t i = 1; i <= n_inters; i++) {
    mask_t mask = mask_ini;
    
    vector_v nom_v ((ray_point - plane_point) * plane_normal);
    vector_v denom_v (ray_direction * plane_normal);
    
    nom_sum   += nom_v.sum(mask);
    denom_sum += denom_v.sum(mask);

    ray_dirs.push_back(ray_direction);
    ray_points.push_back(ray_point);
    result_rds.push_back(vector_v::Zero());

    // n_vec: number of computation vectors that form a geometric vector
    if (i % n_vec == 0) {
      coeff.push_back(nom_sum / denom_sum);

      std::for_each(ray_dirs.begin(), ray_dirs.end(), [&](vector_v &ray_dir){ray_dir.setZeroInverted(mask);});
      for(size_t k = 0; k < result_rds.size(); k++) result_rds[k] += coeff[0]*ray_dirs[k];
      // If you have to sum over subvectors in this case ray_dirs has only one element by construction
      for (size_t j = 0; j < n_sub_loop; j++) {
        // sum over next subvector
        mask = mask.shifted(-4);
        nom_sum   = nom_v.sum(mask);
        denom_sum = denom_v.sum(mask);
        coeff.push_back(nom_sum / denom_sum);

        // restore original data
        ray_dirs[0] = +ray_direction;
        ray_dirs[0].setZeroInverted(mask);
        result_rds[0] += coeff.back()*ray_dirs[0];
      }
      for(size_t k = 0; k < result_rds.size(); k++) result_rds[k] -= ray_points[k];
      results.push_back({.path = std::move(result_rds)
                         .dist = std::move(coeff)});
      // reset for outside loop
      nom_sum = 0; denom_sum = 0;
      ray_dirs.clear();
      ray_points.clear();
      result_rds.clear();
    }

    plane_point.load(pp_ptr + i * stride, Vc::Streaming);
    plane_normal.load(pn_ptr + i * stride, Vc::Streaming);
    // Load different subvectors in raydata too
    if (vector_v::Size/4 < 1) {
      ray_direction.load(rd_ptr + (i % stride) * stride);
      ray_point.load(rp_ptr + (i % stride) * stride);
    }
  }
}*/

/*template<typename scalar_v, typename matrix_t>
void vc_intersect_hybrid(Vector3<scalar_v> &rayVector,
                         Vector3<scalar_v> &rayPoint,
                         aligned::vector<matrix_t> &pns,
                         aligned::vector<matrix_t> &pps,
                         aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results) {

  using scalar_t = typename scalar_v::value_type;

  if (pns.size() != pps.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }

  auto pps_mem = Vc::InterleavedMemoryWrapper<Vector3<scalar_t>, scalar_v>(reinterpret_cast<Vector3<scalar_t>*>(pps.front().data())) ;
  auto pns_mem = Vc::InterleavedMemoryWrapper<Vector3<scalar_t>, scalar_v>(reinterpret_cast<Vector3<scalar_t>*>(pns.front().data())) ;

  for (size_t i = 0; i < pps.size(); i++) {
    scalar_v pps_x, pps_y, pps_z;
    (pps_x, pps_y, pps_z) = pps_mem[i];

    scalar_v pns_x, pns_y, pns_z;
    (pns_x, pns_y, pns_z) = pns_mem[i];

    scalar_v denom_x (rayVector.x * pns_x);
    scalar_v denom_y (rayVector.y * pns_y);
    scalar_v denom_z (rayVector.z * pns_z);

    scalar_v denoms (denom_x + denom_y + denom_z);

    scalar_v nom_x ((rayPoint.x - pps_x) * pns_x);
    scalar_v nom_y ((rayPoint.y - pps_y) * pns_y);
    scalar_v nom_z ((rayPoint.z - pps_z) * pns_z);

    scalar_v coeffs ((nom_x + nom_y + nom_z) / denoms);

    Vector3<scalar_v> path = {.x = (rayPoint.x - coeffs*rayVector.x),
                              .y = (rayPoint.y - coeffs*rayVector.y),
                              .z = (rayPoint.z - coeffs*rayVector.z)};

    results.push_back({.path = std::move(path), .dist = std::move(coeffs)});
  }
}


template<typename scalar_v>
auto vc_intersect_hybrid(Vector3<scalar_v> &rayVector,
                         Vector3<scalar_v> &rayPoint,
                         Vector3<typename scalar_v::value_type> pns,
                         Vector3<typename scalar_v::value_type> pps) {
  using scalar_t = typename scalar_v::value_type;

  if (pns.size() != pps.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }

  auto pps_mem = Vc::InterleavedMemoryWrapper<Vector3<scalar_t>, scalar_v>(&pps);
  auto pns_mem = Vc::InterleavedMemoryWrapper<Vector3<scalar_t>, scalar_v>(&pns);

  scalar_v pps_x, pps_y, pps_z;
  (pps_x, pps_y, pps_z) = pps_mem[0];

  scalar_v pns_x, pns_y, pns_z;
  (pns_x, pns_y, pns_z) = pns_mem[0];

  scalar_v denom_x (rayVector.x * pns_x);
  scalar_v denom_y (rayVector.y * pns_y);
  scalar_v denom_z (rayVector.z * pns_z);

  scalar_v denoms (denom_x + denom_y + denom_z);

  scalar_v nom_x ((rayPoint.x - pps_x) * pns_x);
  scalar_v nom_y ((rayPoint.y - pps_y) * pns_y);
  scalar_v nom_z ((rayPoint.z - pps_z) * pns_z);

  scalar_v coeffs ((nom_x + nom_y + nom_z) / denoms);

  Vector3<scalar_v> path = {.x = (rayPoint.x - coeffs*rayVector.x),
                            .y = (rayPoint.y - coeffs*rayVector.y),
                            .z = (rayPoint.z - coeffs*rayVector.z)};

  intersection<scalar_v, Vector3<scalar_v>> res = {.path = std::move(path), .dist = std::move(coeffs)};
  return res;
}*/


template<typename scalar_v>
void vc_intersect_hybrid(ray_data<Vector3<scalar_v>>& ray,
                         plane_data<aligned::vector<Vector3<typename scalar_v::value_type>>> &planes,
                         aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results) {

  using scalar_t = typename scalar_v::value_type;

  /*if (pns_struct.size() != pps_struct.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }*/

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

  using scalar_t = typename scalar_v::value_type;

  /*if (pns_struct.size() != pps_struct.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }*/

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

  intersection<scalar_v, Vector3<scalar_v>> res = {.path = std::move(path), .dist = std::move(coeffs)};
  return res;
}


template<typename scalar_v, typename vector_s, size_t kDIM = 3>
void vc_intersect_horiz(ray_data<Vector3<scalar_v>>& ray,
                        plane_data<aligned::vector<MatrixV<vector_s>>>& planes,
                        aligned::vector<intersection<scalar_v, Vector3<scalar_v>>> &results,
                        size_t padding = 0) {
  //TODO: make constexpr
  /*if (planePoints.size() != planeNormals.size()) {
    std::cerr << "Error: Different size of input collections (plane points and normals)" << std::endl;
    return;
  }*/

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

    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction.x, -ray.point.x), 
                              .y = Vc::fma(coeffs, ray.direction.y, -ray.point.y), 
                              .z = Vc::fma(coeffs, ray.direction.z, -ray.point.z)};

    results.emplace_back(intersection<scalar_v, Vector3<scalar_v>>{.path = std::move(path), .dist = std::move(coeffs)});
  }
}


template<typename scalar_v, typename data_ptr_t, size_t kDIM = 3>
auto vc_intersect_horiz(ray_data<Vector3<scalar_v>>& ray,
                        data_ptr_t pl_normals_ptr,
                        data_ptr_t pl_points_ptr,
                        size_t& offset) {
    /*if (pl_normals_ptr == nullptr || pl_points_ptr == nullptr) {
      std::cerr << "Passed invalid data collection pointer to intersection" << std::endl;
      return intersection<scalar_v, Vector3<scalar_v>>{};
    }*/

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

    Vector3<scalar_v> path = {.x = Vc::fma(coeffs, ray.direction.x, -ray.point.x), 
                              .y = Vc::fma(coeffs, ray.direction.y, -ray.point.y), 
                              .z = Vc::fma(coeffs, ray.direction.z, -ray.point.z)};

    intersection<scalar_v, Vector3<scalar_v>> res = {.path = std::move(path), .dist = std::move(coeffs)};
    return std::move(res);
}


} // namespace vec_intr