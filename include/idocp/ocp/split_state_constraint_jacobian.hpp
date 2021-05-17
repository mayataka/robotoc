#ifndef IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HPP_ 
#define IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

///
/// @class SplitStateConstraintJacobian
/// @brief Jacobian of SwitchingConstraint. 
///
class SplitStateConstraintJacobian {
public:
  ///
  /// @brief Constructs a Jacobian of the switching constraint.
  /// @param[in] robot Robot model. 
  ///
  SplitStateConstraintJacobian(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitStateConstraintJacobian();

  ///
  /// @brief Destructor. 
  ///
  ~SplitStateConstraintJacobian();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitStateConstraintJacobian(const SplitStateConstraintJacobian&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitStateConstraintJacobian& operator=(
      const SplitStateConstraintJacobian&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitStateConstraintJacobian(
      SplitStateConstraintJacobian&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitStateConstraintJacobian& operator=(
      SplitStateConstraintJacobian&&) noexcept = default;

  ///
  /// @brief Set the dimension of the impulse. 
  /// @param[in] impulse_status Impulse status. 
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set the dimension of the impulse 0. 
  ///
  void setImpulseStatus();

  int dimi() const;

  Eigen::Block<Eigen::MatrixXd> Phia();

  const Eigen::Block<const Eigen::MatrixXd> Phia() const;

  Eigen::Block<Eigen::MatrixXd> Phix();

  const Eigen::Block<const Eigen::MatrixXd> Phix() const;

  Eigen::Block<Eigen::MatrixXd> Phiq();

  const Eigen::Block<const Eigen::MatrixXd> Phiq() const;

  Eigen::Block<Eigen::MatrixXd> Phiv();

  const Eigen::Block<const Eigen::MatrixXd> Phiv() const;

  Eigen::Block<Eigen::MatrixXd> Phiu();

  const Eigen::Block<const Eigen::MatrixXd> Phiu() const;

  Eigen::MatrixXd dintegrate_dq;

  Eigen::MatrixXd dintegrate_dv;

  bool isApprox(const SplitStateConstraintJacobian& other) const;

  bool hasNaN() const;

private:
  Eigen::MatrixXd Phia_full_, Phix_full_, Phiu_full_;
  int dimv_, dimx_, dimu_, dimi_;

};

} // namespace idocp

#include "idocp/ocp/split_state_constraint_jacobian.hxx"

#endif // IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HPP_