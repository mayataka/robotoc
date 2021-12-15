#ifndef ROBOTOC_DISCRETIZATION_METHOD_HPP_ 
#define ROBOTOC_DISCRETIZATION_METHOD_HPP_

namespace robotoc {

/// 
/// @enum DiscretizationMethod
/// @brief Discretization method of the hybrid optimal control problem.
///
enum class DiscretizationMethod {
  GridBased,
  PhaseBased
};

} // namespace robotoc

#endif // ROBOTOC_DISCRETIZATION_METHOD_HPP_ 