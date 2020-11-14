#ifndef IDOCP_HYBRID_CONTAINER_HPP_
#define IDOCP_HYBRID_CONTAINER_HPP_

#include <vector>

namespace idocp {

///
/// @class hybrid_container
/// @brief A container that is useful to formulate the hybrid optimal control 
/// problem. This container has the standard data (with Type), data for lift 
/// stages (with Type), and data for impulse stages (with ImpulseType).
/// @tparam Type The type name of the standard data type.
/// @tparam ImpulseType The type name of the impulse data type.
/// 
///
template <typename Type, typename ImpulseType>
struct hybrid_container {

  hybrid_container(const int N, const int max_num_impulse_stages) 
    : data(N, Type()), 
      impulse(max_num_impulse_stages, ImpulseType()),
      lift(max_num_impulse_stages, Type()) {
  }

  hybrid_container(const int N, const Type& obj, 
                   const int max_num_impulse_stages, 
                   const ImpulseType& impulse_obj) 
    : data(N, obj), 
      impulse(max_num_impulse_stages, impulse_obj),
      lift(max_num_impulse_stages, obj) {
  }

  hybrid_container(const int N) 
    : data(N, Type()), 
      impulse(),
      lift() {
  }

  hybrid_container(const int N, const Type& obj) 
    : data(N, obj), 
      impulse(),
      lift() {
  }

  hybrid_container() 
    : data(), 
      impulse(),
      lift() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  hybrid_container(const hybrid_container&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  hybrid_container& operator=(const hybrid_container&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  hybrid_container(hybrid_container&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  hybrid_container& operator=(hybrid_container&&) noexcept = default;

  ///
  /// @brief Overload operator[] to access the standard data as std::vector. 
  ///
  Type& operator[] (const int i) {
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const Type& operator[] (const int i) const {
    return data[i];
  }

  std::vector<Type> data, lift;
  std::vector<ImpulseType> impulse;
};
  
} // namespace idocp

#endif // IDOCP_HYBRID_CONTAINER_HPP_