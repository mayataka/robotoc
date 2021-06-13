#ifndef IDOCP_SPLIT_SWITCHING_CONSTRAINT_RESIDUAL_HPP_ 
#define IDOCP_SPLIT_SWITCHING_CONSTRAINT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class SplitSwitchingConstraintResidual
/// @brief KKT residual w.r.t. the switching constraint split into each time stage. 
///
class SplitSwitchingConstraintResidual {
public:
  ///
  /// @brief Construct a split KKT residual.
  /// @param[in] robot Robot model. 
  ///
  SplitSwitchingConstraintResidual(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitSwitchingConstraintResidual();

  ///
  /// @brief Destructor. 
  ///
  ~SplitSwitchingConstraintResidual();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitSwitchingConstraintResidual(
      const SplitSwitchingConstraintResidual&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitSwitchingConstraintResidual& operator=(
      const SplitSwitchingConstraintResidual&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitSwitchingConstraintResidual(
      SplitSwitchingConstraintResidual&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitSwitchingConstraintResidual& operator=(
      SplitSwitchingConstraintResidual&&) noexcept = default;

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse, to zero.
  ///
  void setImpulseStatus();

  ///
  /// @brief Residual in the switching constraint.
  /// @return Reference to the residual in the switching constraints. 
  /// Size is ImpulseStatus::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> P();

  ///
  /// @brief const version of SplitSwitchingConstraintResidual::P().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> P() const;

  ///
  /// @brief Predicted configuration used for computing the switching 
  /// constraint. Size is Robot::dimq().
  ///
  Eigen::VectorXd q;

  ///
  /// @brief Predicted difference of the configuration used for computing the 
  /// switching constraint. Size is Robot::dimv().
  ///
  Eigen::VectorXd dq;

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual of the switching constraint. 
  /// @return Squared norm of the KKT residual in the switching constraint.
  ///
  double squaredNormKKTResidual() const;

  ///
  /// @brief Returns l1-norm of the constraint violation, that is, the primal
  /// residual in the switchign constraint. 
  /// @return l1-norm of the constraint violation.
  ///
  double l1NormConstraintViolation() const;

  ///
  /// @brief Sets the split KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of the stack of impulse forces.
  ///
  int dimi() const;

  ///
  /// @brief Checks dimensional consistency of each component. 
  /// @return true if the dimension is consistent. false if not.
  ///
  bool isDimensionConsistent() const;

  ///
  /// @brief Checks the equivalence of two SplitSwitchingConstraintResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitSwitchingConstraintResidual& other) const;

  ///
  /// @brief Checks this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

private:
  Eigen::VectorXd P_full_;
  int dimq_, dimv_, dimi_;

};

} // namespace idocp 

#include "idocp/ocp/split_switching_constraint_residual.hxx"

#endif // IDOCP_SPLIT_SWITCHING_CONSTRAINT_RESIDUAL_HPP_ 