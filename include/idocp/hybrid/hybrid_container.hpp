#ifndef IDOCP_HYBRID_CONTAINER_HPP_
#define IDOCP_HYBRID_CONTAINER_HPP_

#include <vector>

namespace idocp {

///
/// @class hybrid_container
/// @brief A container that is useful to formulate the hybrid optimal control 
/// problem. This container has the standard data (with Type), data for lift 
/// stages (with Type), data for aux stages (with Type), and data for impulse 
/// stages (with ImpulseType).
/// @tparam Type The type name of the standard data type.
/// @tparam ImpulseType The type name of the impulse data type.
/// 
///
template <typename Type, typename ImpulseType>
struct hybrid_container {
  ///
  /// @brief Construct the standard data, impulse data, and lift data. 
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  ///
  hybrid_container(const int N, const int N_impulse) 
    : data(N, Type()), 
      aux(N_impulse, Type()), 
      lift(N_impulse, Type()),
      impulse(N_impulse, ImpulseType()) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  /// @param[in] obj An object of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  /// @param[in] impulse_obj An object of the impulse data.
  ///
  hybrid_container(const int N, const Type& obj, const int N_impulse, 
                   const ImpulseType& impulse_obj) 
    : data(N, obj), 
      aux(N_impulse, obj), 
      lift(N_impulse, obj),
      impulse(N_impulse, impulse_obj) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  ///
  hybrid_container(const int N) 
    : data(N, Type()), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  /// @param[in] obj An object of the standard data.
  ///
  hybrid_container(const int N, const Type& obj) 
    : data(N, obj), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Default Constructor.
  ///
  hybrid_container() 
    : data(), 
      aux(),
      lift(),
      impulse() {
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
  /// @brief Overload operator[] to access the standard data, i.e., 
  /// hybrid_container::data as std::vector. 
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

  std::vector<Type> data, aux, lift;
  std::vector<ImpulseType> impulse;
};
  
} // namespace idocp

#endif // IDOCP_HYBRID_CONTAINER_HPP_