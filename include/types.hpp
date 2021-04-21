#include <limits>
#include <type_traits>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

#include <Vc/Vc>

//#include <vecmem/containers/vector.hpp>
//#include <vecmem/memory/host_memory_resource.hpp>

#pragma once

namespace vec_intr {

#ifdef PRECISION_DOUBLE
using Scalar = double;
using Scalar_v = Vc::double_v;
#else
using Scalar = float;
using Scalar_v = Vc::float_v;
#endif

using Index_v  = Scalar_v::IndexType;

using Vector3_s  = Eigen::Matrix<Scalar, 3, 1>;
using Vector4_s  = Eigen::Matrix<Scalar, 4, 1>;
using VectorV_s  = Eigen::Matrix<Scalar, Scalar_v::Size, 1>;
using Transform4 = Eigen::Transform<Scalar, 4, Eigen::Affine>;


// prevent false cache sharing
#ifdef __cpp_lib_hardware_interference_size
    using std::hardware_constructive_interference_size;
    using std::hardware_destructive_interference_size;
#else
    // Lucky guess │ __cacheline_aligned │ L1_CACHE_BYTES │ L1_CACHE_SHIFT │ ...
    constexpr std::size_t hardware_constructive_interference_size
        = 2 * sizeof(std::max_align_t);
    constexpr std::size_t hardware_destructive_interference_size
        = 2 * sizeof(std::max_align_t);
#endif

// Check for max alignment!!e.g. Vc::VectorAlignment 
//constexpr size_t alignment = hardware_constructive_interference_size;
constexpr size_t alignment = alignof(Scalar_v);
namespace aligned {
  // Allow Vc to get alignment right
  //template <typename T, typename Allocator = std::allocator<T>>
  //template <typename T, typename Allocator = Eigen::aligned_allocator<T>>
  template <typename T, typename Allocator = Vc::Allocator<T>>
  // Add subscript operator to allow for gather operations from AoS
  using vector = Vc::Common::AdaptSubscriptOperator<std::vector<T, Allocator>>;

  template<size_t kDIM>
  using mem_t = Vc::Memory<Scalar_v, kDIM>;

  template<typename vector_s, typename mem_v>
  using storage_t = std::aligned_union_t<1, vector_s, mem_v>;

  template<typename vector_t, typename memory_t>
  union data_container {vector_t vector;
                        memory_t memory;};
} //namespace aligned

// Use in AoS
template<typename data_t>
struct Vector3
{
  data_t x, y, z;
};

template<typename data_t>
struct Vector4
{
  data_t x, y, z, t;
};

// Use in AoS for vertical vectorization
// Keep the geometrical vectors as Vc vectors (vertical vect.)
// Keep the plane points extra to make the struct alignment easier
template<typename data_t>
struct alignas(alignment) ray_data
{
  data_t direction, point;
};

// Use in AoS for vertical vectorization
// Keep the geometrical vectors as Vc vectors (vertical vect.)
// Keep the plane points extra to make the struct alignment easier
template<typename data_t>
struct alignas(alignment) plane_data
{
  data_t normals, points;
};

// detray style intersection structure to be filled as result type
template <typename scalar_t, typename vector_t = Vector3_s>
struct alignas(alignment) intersection {
  vector_t path;
  scalar_t dist = std::numeric_limits<scalar_t>::infinity();
};

// define type that holds input data: has to be memory layout compatible with the Vc Vector type, i.e. 
//                                    has to have same ABI as array of T with rows*cols many entries
// To be specialized by the LA plugins
/*template<typename T, class Enable = void>
struct data_trait {
  using type = void;
  using value_type = void;

  static constexpr bool is_vec_layout = false;
  static constexpr bool is_std_layout = false;
  static constexpr bool is_same_type  = false;
  static constexpr bool is_vec_dim    = false;
};*/

//---------------------------------------------------Plugin specific

// define type that holds input data: has to be memory layout compatible with the Vc Vector type, i.e. 
//                                    has to have same ABI as array of T with rows*cols many entries
// To be specialized by the LA plugins
template<typename Derived>
struct alignas(alignment) MatrixV {
  using value_type = typename Eigen::DenseBase<Derived>::Scalar;
  using obj_type   = Derived;
  obj_type obj;
  //void init() {obj = decltype(obj)::Random();}
  constexpr size_t n_elemts() {return obj.rows() * obj.cols();}
  const value_type* data() {return obj.data();}
  constexpr size_t padding() {return (alignment - n_elemts() * sizeof(value_type)) / sizeof(value_type);}
};

// Affine transform is not derived from Eigen::DenseBase
template<>
struct alignas(alignment) MatrixV<Transform4> {
  using value_type = typename Transform4::Scalar;
  using obj_type   = Transform4;
  Transform4 obj;
  const size_t n_elemts() {return obj.rows() * obj.cols();}
  const value_type* data() {return obj.data();}
  const size_t padding() {return (alignment - n_elemts() * sizeof(value_type)) / sizeof(value_type);}
};


template<typename Derived>
struct data_trait {
  using type  = MatrixV<Derived>;
  using value_type = typename type::value_type;
  static const size_t size = sizeof(type)/sizeof(value_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type>;
  static constexpr bool is_same_type  = std::is_same_v<value_type, Scalar_v::value_type>;
  static constexpr bool is_vec_dim    = !(Scalar_v::Size > size ? Scalar_v::Size % size : size % Scalar_v::Size);
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};

// Affine transform is not derived from Eigen::DenseBase
template<>
struct data_trait<Transform4> {
  using type  = MatrixV<Transform4>;
  using value_type = Transform4::Scalar;
  static const size_t size = sizeof(type)/sizeof(value_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type> 
                                        && std::is_standard_layout_v<Transform4::VectorType>;
  static constexpr bool is_same_type  = std::is_same_v<value_type, Scalar_v::value_type>;
  static constexpr bool is_vec_dim    = Scalar_v::Size > size ? Scalar_v::Size % size == 0 : size % Scalar_v::Size == 0;
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};

} //namespace vec_intr