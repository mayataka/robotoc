#ifndef IDOCP_HYBRID_CONTAINER_HPP_
#define IDOCP_HYBRID_CONTAINER_HPP_

#include <vector>
#include <cassert>

#include "idocp/robot/robot.hpp"


namespace idocp {

namespace internal {
///
/// @class EmptyType
/// @brief Only used for default template parameter of hybrid_container.
/// This class does not do anything.
///
class EmptyType {
public:
  EmptyType(const Robot&) {}
  EmptyType() {}
  ~EmptyType() {}
};
} // namespace internal


///
/// @class hybrid_container
/// @brief A container used for formulating the hybrid optimal control problem. 
/// This container has data for time stages (with Type), data for lift stages 
/// (with Type), data for aux stages (with Type), data for impulse stages
/// (with ImpulseType), and data for switching ingredients.
/// @tparam Type The type name of the standard data type.
/// @tparam ImpulseType The type name of the impulse data type. Defalt is 
/// internal::EmptyType.
/// @tparam SwitchingType The type name of the switching data type. Defalt is 
/// internal::EmptyType.
///
template <typename Type, typename ImpulseType=internal::EmptyType, 
          typename SwitchingType=internal::EmptyType>
class hybrid_container {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data. Default is 0.
  ///
  hybrid_container(const Robot& robot, const int N, const int N_impulse=0) 
    : data(N+1, Type(robot)), 
      aux(N_impulse, Type(robot)), 
      lift(N_impulse, Type(robot)),
      impulse(N_impulse, ImpulseType(robot)),
      switching(N_impulse, SwitchingType(robot)) {
  }

  ///
  /// @brief Default Constructor.
  ///
  hybrid_container() 
    : data(), 
      aux(),
      lift(),
      impulse(),
      switching() {
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
  /// @brief Overload operator[] to access the hybrid_container::data as 
  /// std::vector. 
  ///
  Type& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const Type& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<Type> data, aux, lift;
  std::vector<ImpulseType> impulse;
  std::vector<SwitchingType> switching;
};

} // namespace idocp

#endif // IDOCP_HYBRID_CONTAINER_HPP_