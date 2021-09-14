#ifndef IDOCP_SPLIT_SWITCHING_CONSTRAINT_JACOBIAN_HPP_ 
#define IDOCP_SPLIT_SWITCHING_CONSTRAINT_JACOBIAN_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class SplitSwitchingConstraintJacobian
/// @brief The KKT matrix split into a time stage.
///
class SplitSwitchingConstraintJacobian {
public:
  ///
  /// @brief Construct a split KKT matrix.
  /// @param[in] robot Robot model. 
  ///
  SplitSwitchingConstraintJacobian(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitSwitchingConstraintJacobian();

  ///
  /// @brief Destructor. 
  ///
  ~SplitSwitchingConstraintJacobian();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitSwitchingConstraintJacobian(
      const SplitSwitchingConstraintJacobian&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitSwitchingConstraintJacobian& operator=(
      const SplitSwitchingConstraintJacobian&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitSwitchingConstraintJacobian(
      SplitSwitchingConstraintJacobian&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitSwitchingConstraintJacobian& operator=(
      SplitSwitchingConstraintJacobian&&) noexcept = default;

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulses.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse, to zero.
  ///
  void setImpulseStatus();

  ///
  /// @brief Jacobian of the original contact position constraint w.r.t. q. 
  /// @return Reference to the Jacobian. 
  /// Size is ImpulseStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Pq();

  ///
  /// @brief const version of SplitSwitchingConstraintJacobian::Pq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Pq() const;

  ///
  /// @brief Jacobian of the swithcing constraint w.r.t. x. 
  /// @return Reference to the Jacobian. 
  /// Size is ImpulseStatus::dimf() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Phix();

  ///
  /// @brief const version of SplitSwitchingConstraintJacobian::Phix().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Phix() const;

  ///
  /// @brief Jacobian of the swithcing constraint w.r.t. q. 
  /// @return Reference to the Jacobian. 
  /// Size is ImpulseStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Phiq();

  ///
  /// @brief const version of SplitSwitchingConstraintJacobian::Phiq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Phiq() const;

  ///
  /// @brief Jacobian of the swithcing constraint w.r.t. v. 
  /// @return Reference to the Jacobian. 
  /// Size is ImpulseStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Phiv();

  ///
  /// @brief const version of SplitSwitchingConstraintJacobian::Phiv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Phiv() const;

  ///
  /// @brief Jacobian of the swithcing constraint w.r.t. a. 
  /// @return Reference to the Jacobian. 
  /// Size is ImpulseStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Phia();

  ///
  /// @brief const version of SplitSwitchingConstraintJacobian::Phia().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Phia() const;

  ///
  /// @brief Jacobian of the swithcing constraint w.r.t. u. 
  /// @return Reference to the Jacobian. 
  /// Size is ImpulseStatus::dimf() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Phiu();

  ///
  /// @brief const version of SplitSwitchingConstraintJacobian::Phiu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Phiu() const;

  ///
  /// @brief Jacobian of the swithcing constraint w.r.t. the switching time. 
  /// @return Reference to the time Jacobian vector. 
  /// Size is ImpulseStatus::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Phit();

  ///
  /// @brief const version of SplitSwitchingConstraintJacobian::Phit().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Phit() const;

  ///
  /// @brief Set the all components zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of the stack of impulse forces.
  ///
  int dimi() const;

  ///
  /// @brief Checks the equivalence of two SplitSwitchingConstraintJacobian.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitSwitchingConstraintJacobian& other) const;

  ///
  /// @brief Checks this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

private:
  Eigen::MatrixXd Pq_full_, Phix_full_, Phia_full_, Phiu_full_;
  Eigen::VectorXd Phit_full_;
  bool has_floating_base_;
  int dimv_, dimx_, dimu_, dimi_;

};

} // namespace idocp 

#include "idocp/ocp/split_switching_constraint_jacobian.hxx"

#endif // IDOCP_SPLIT_SWITCHING_CONSTRAINT_JACOBIAN_HPP_ 