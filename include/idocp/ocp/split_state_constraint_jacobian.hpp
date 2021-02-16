#ifndef IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HPP_ 
#define IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

class SplitStateConstraintJacobian {
public:
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

  Eigen::VectorXd q;

  Eigen::VectorXd dq;

  Eigen::MatrixXd dintegrate_dq;

  Eigen::MatrixXd dintegrate_dv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd Phia_full_, Phix_full_, Phiu_full_;
  int dimv_, dimx_, dimu_, dimi_;

};

} // namespace idocp

#include "idocp/ocp/split_state_constraint_jacobian.hxx"

#endif // IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HPP_