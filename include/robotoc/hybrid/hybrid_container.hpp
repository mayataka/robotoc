#ifndef ROBOTOC_HYBRID_CONTAINER_HPP_
#define ROBOTOC_HYBRID_CONTAINER_HPP_

#include <vector>
#include <cassert>
#include <iostream>

#include "robotoc/robot/robot.hpp"


namespace robotoc {

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

  ///
  /// @brief Does not print anything.
  ///
  friend std::ostream& operator<<(std::ostream& os, const EmptyType& obj) {
    return os;
  }
};
} // namespace internal


///
/// @class hybrid_container
/// @brief A container used for formulating the hybrid optimal control problem. 
/// This container has data for time stages with Type, data for lift stages 
/// (additional time stages just after the lift events) (with Type), data for 
/// the auxiliary stages (additional time stages just after the impulse events)
/// (with Type), data for impulse stages (additional time stages at the impulse 
/// events) (with ImpulseType), and data for switching constraints (with 
/// SwitchingType).
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
  /// @param[in] N Number of the discretization grids of the horizon except for 
  /// the discrete events. Must be positive.
  /// @param[in] max_num_each_discrete_events Maximum possible number of the 
  /// each discrete events on the horizon. OCP data for impulse and lift events 
  /// are constructed according to this value. Default is 0. 
  /// Must be non-negative.
  ///
  hybrid_container(const Robot& robot, const int N, 
                   const int max_num_each_discrete_events=0) 
    : data(N+1, Type(robot)), 
      aux(max_num_each_discrete_events, Type(robot)), 
      lift(max_num_each_discrete_events, Type(robot)),
      impulse(max_num_each_discrete_events, ImpulseType(robot)),
      switching(max_num_each_discrete_events, SwitchingType(robot)) {
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
  /// @brief Resize the data without reallocating all the data. 
  /// @param[in] robot Robot model.
  /// @param[in] max_num_each_discrete_events Maximum possible number of the 
  /// each discrete events on the horizon. OCP data for impulse and lift events 
  /// are constructed according to this value. Must be non-negative.
  /// @remark robot model must be the same as one in the argument of the 
  /// constructor.
  ///
  void resize(const Robot& robot, const int max_num_each_discrete_events) {
    assert(max_num_each_discrete_events >= 0);
    assert(impulse.size() == aux.size());
    assert(impulse.size() == switching.size());
    while (max_num_each_discrete_events > impulse.size()) {
      impulse.emplace_back(robot);
      aux.emplace_back(robot);
      switching.emplace_back(robot);
    }
    while (max_num_each_discrete_events > lift.size()) {
      lift.emplace_back(robot);
    }
  }

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

  ///
  /// @brief Data for the time stages.
  ///
  std::vector<Type> data;

  ///
  /// @brief Data for the auxiliary stages 
  /// (additional time stages just after the impulse events).
  ///
  std::vector<Type> aux;

  ///
  /// @brief Data for the lift stages 
  /// (additional time stages just after the lift events).
  ///
  std::vector<Type> lift;

  ///
  /// @brief Data for the impulse stages 
  /// (additional time stages at the impulse events).
  ///
  std::vector<ImpulseType> impulse;

  ///
  /// @brief Data for the switching constraints. 
  ///
  std::vector<SwitchingType> switching;

  ///
  /// @brief Displays the elements of the container onto a ostream.
  ///
  void disp(std::ostream& os) const {
    const int N = data.size() - 1;
    for (int i=0; i<=N; ++i) {
      os << "  stage: " << i << std::endl;
      os << data[i] << std::endl;
    }
    const int max_num_impulses = impulse.size();
    for (int i=0; i<max_num_impulses; ++i) {
      os << "  impulse: " << i << std::endl;
      os << impulse[i] << std::endl;
    }
    for (int i=0; i<max_num_impulses; ++i) {
      os << "  aux: " << i << std::endl;
      os << aux[i] << std::endl;
    }
    for (int i=0; i<max_num_impulses; ++i) {
      os << "  switching: " << i << std::endl;
      os << switching[i] << std::endl;
    }
    const int max_num_lifts = lift.size();
    for (int i=0; i<max_num_lifts; ++i) {
      os << "  lift: " << i << std::endl;
      os << lift[i] << std::endl;
    }
  }

};

} // namespace robotoc

#endif // ROBOTOC_HYBRID_CONTAINER_HPP_