#ifndef IDOCP_FLOATING_BASE_HPP_
#define IDOCP_FLOATING_BASE_HPP_

#include <vector>

#include "Eigen/Core"
#include "pinocchio/multibody/model.hpp"


namespace idocp {

///
/// @class FloatingBase
/// @brief Floating base model.
///
class FloatingBase {
public:
  ///
  /// @brief Construct floating base model.
  /// @param[in] model The pinocchio model. Before call this constructor, 
  /// pinocchio model must be initialized, e.g., by pinocchio::buildModel().
  ///
  FloatingBase(const pinocchio::Model& model);

  ///
  /// @brief Default constructor. 
  ///
  FloatingBase();

  ///
  /// @brief Destructor. 
  ///
  ~FloatingBase();

  ///
  /// @brief Default copy constructor. 
  ///
  FloatingBase(const FloatingBase&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  FloatingBase& operator=(const FloatingBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  FloatingBase(FloatingBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  FloatingBase& operator=(FloatingBase&&) noexcept = default;

  ///
  /// @brief Returns the dimensiton of the generalized torques corresponding to 
  /// the passive joints.
  /// @return The dimensiton of the generalized torques corresponding to the 
  //// passive joints.
  /// 
  int dim_passive() const;

  ///
  /// @brief Returns joint indices of the passive joints.
  /// @return Joitn indices of the passive joints.
  /// 
  std::vector<int> passive_joint_indices() const;

  ///
  /// @brief Returns true if the robot has a floating base and false if not.
  /// @returns true if the robot has a floating base and false if not.
  /// 
  bool has_floating_base() const;

private:
  bool has_floating_base_;
  std::vector<int> passive_joint_indices_;
  int dimv_, dim_passive_;
};

} // namespace idocp

#include "idocp/robot/floating_base.hxx"

#endif // IDOCP_FLOATING_BASE_HPP_