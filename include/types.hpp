 /**
 * author: joana.niermann@cern.ch
 **/
#include <limits>
#include <type_traits>

#include <Eigen/Dense>
#include <Eigen/Geometry>

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


template<typename data_t>
struct Vector3;

using Vector3_s    = Eigen::Matrix<Scalar, 3, 1>;
using Vector4_s    = Eigen::Matrix<Scalar, 4, 1>;
using VectorV_s    = Eigen::Matrix<Scalar, Scalar_v::Size, 1>;
using Vector3_v    = Vector3<Scalar_v>;
using Transform4_s = Eigen::Transform<Scalar, 4, Eigen::Affine>;

// Check for max alignment!!e.g. Vc::VectorAlignment 
//constexpr size_t alignment = std::hardware_constructive_interference_size;
constexpr size_t alignment = alignof(Scalar_v);
namespace aligned {
  // Allow Vc to get alignment right
  //template <typename T, typename Allocator = std::allocator<T> >
  //template <typename T, typename Allocator = Eigen::aligned_allocator<T> >
  template <typename T, typename Allocator = Vc::Allocator<T> >
  // Add subscript operator to allow for gather operations from AoS
  using vector = Vc::Common::AdaptSubscriptOperator<std::vector<T, Allocator> >;

  template<size_t kDIM>
  using mem_t = Vc::Memory<Scalar_v, kDIM>;

  // TODO use placement new rather
  template<typename vector_s, typename mem_v>
  using storage_t = std::aligned_union_t<1, vector_s, mem_v>;

  template<typename vector_t, typename memory_t>
  union data_container {vector_t vector;
                        memory_t memory;};
} //namespace aligned

// Needed in AoS
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

// Assume 4x4 transformation matrix (placement) as input
template<typename data_t>
struct Transform4
{
  // the rest of the 4x4 matrix (avoid padding)
  Vector4<data_t> vec0, vec1;
  // plane normal
  Vector4<data_t> normal;
  // plane translation
  Vector4<data_t> translation;
};

// Assume 4x4 transformation (placement) as input 
// Eigen compact format (3x3 in memory)
template<typename data_t>
struct Transform3
{
  // the rest of the 3x3 matrix (avoid padding)
  Vector3<data_t> vec0, vec1;
  // plane normal
  Vector3<data_t> normal;
  // plane translation
  Vector3<data_t> translation;
};

// Convenience types
template<typename data_t>
struct alignas(alignment) ray_data
{
  data_t direction, point;
};

template<typename data_t>
struct alignas(alignment) plane_data
{
  data_t normals, points;
};

// detray style intersection structure to be filled as result type
template <typename scalar_t, typename vector_t = Vector3_s>
struct /*alignas(alignment)*/ intersection {
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

// Small-time wrapper for a uniform interface, will be taken care of by the LA plugin
// later
template<typename Derived>
struct alignas(alignment) MatrixV {
  using scalar_type = typename Eigen::DenseBase<Derived>::Scalar;
  using vec_type = typename Eigen::DenseBase<Derived>::Scalar;
  using type   = Derived;

  type obj;

  type& operator()() {return obj;}

  //void init() {obj = decltype(obj)::Random();}
  constexpr size_t n_elemts() {return obj.rows() * obj.cols();}
  const scalar_type* data() {return obj.data();}
  constexpr size_t padding() {return (alignment - (n_elemts() * sizeof(scalar_type)) % alignment) % alignment / sizeof(scalar_type);}
};

// Affine transform is not derived from Eigen::DenseBase
template<>
struct alignas(alignment) MatrixV<Transform4_s> {
  using scalar_type = typename Transform4_s::Scalar;
  using vec_type = typename Transform4_s::Scalar;
  using type   = Transform4_s;

  Transform4_s obj;

  type& operator()() {return obj;}
   
  size_t n_elemts() {return obj.rows() * obj.cols();}
  const scalar_type* data() {return obj.data();}
  size_t padding() {return (alignment - (n_elemts() * sizeof(scalar_type)) % alignment) % alignment / sizeof(scalar_type);}
};

// 3D Structure of vector data
template<>
struct alignas(alignment) MatrixV<Vector3<Scalar_v> > {
  using scalar_type = Scalar_v::value_type;
  using vec_type = Scalar_v;
  using type   = Vector3<Scalar_v>;

  Vector3<Scalar_v> obj;

  //MatrixV() {}
  type& operator()() {return obj;}

  size_t n_elemts() {return 3*Scalar_v::Size;}
  const scalar_type* data() {return reinterpret_cast<const scalar_type*>(&obj.x);}
  size_t padding() {return (alignment - (n_elemts() * sizeof(scalar_type)) % alignment) % alignment / sizeof(scalar_type);}
};

// Eigen specific types
template<typename Derived>
struct data_trait {
  using type  = MatrixV<Derived>;
  using scalar_type = typename type::scalar_type;
  using vec_type = typename type::vec_type;
  static const size_t size = sizeof(type)/sizeof(scalar_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type>;
  static constexpr bool is_same_type  = std::is_same_v<scalar_type, Scalar_v::value_type>;
  static constexpr bool is_vec_dim    = !(Scalar_v::Size > size ? Scalar_v::Size % size : size % Scalar_v::Size);
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};

// Affine transform is not derived from Eigen::DenseBase
template<>
struct data_trait<Transform4_s> {
  using type  = MatrixV<Transform4_s>;
  using scalar_type = Transform4_s::Scalar;
  using vec_type = typename type::vec_type;
  static const size_t size = sizeof(type)/sizeof(scalar_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type> 
                                        && std::is_standard_layout_v<Transform4_s::VectorType>;
  static constexpr bool is_same_type  = std::is_same_v<scalar_type, Scalar_v::value_type>;
  static constexpr bool is_vec_dim    = Scalar_v::Size > size ? Scalar_v::Size % size == 0 : size % Scalar_v::Size == 0;
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};

template<>
struct data_trait<Vector3<Scalar_v> > {
  using type  = MatrixV<Vector3<Scalar_v> >;
  using scalar_type = Scalar_v::value_type;
  using vec_type = typename type::vec_type;
  static const size_t size = sizeof(type)/sizeof(scalar_type);

  static constexpr bool is_std_layout = std::is_standard_layout_v<type>;
  static constexpr bool is_same_type  = std::is_same_v<scalar_type, Scalar_v::value_type>;
  static constexpr bool is_vec_dim    = Scalar_v::Size > size ? Scalar_v::Size % size == 0 : size % Scalar_v::Size == 0;
  static constexpr bool is_vec_layout = is_std_layout && is_same_type && is_vec_dim;
};

} //namespace vec_intr