#ifndef ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HPP_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"


namespace robotoc {

///
/// @class SwitchingConstraintResidual
/// @brief KKT residual w.r.t. the switching constraint split into each time stage. 
///
class SwitchingConstraintResidual {
public:
  ///
  /// @brief Construct a split KKT residual.
  /// @param[in] robot Robot model. 
  ///
  SwitchingConstraintResidual(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SwitchingConstraintResidual();

  ///
  /// @brief Destructor. 
  ///
  ~SwitchingConstraintResidual();

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingConstraintResidual(const SwitchingConstraintResidual&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SwitchingConstraintResidual& operator=(
      const SwitchingConstraintResidual&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingConstraintResidual(SwitchingConstraintResidual&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingConstraintResidual& operator=(
      SwitchingConstraintResidual&&) noexcept = default;

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
  /// @brief const version of SwitchingConstraintResidual::P().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> P() const;

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual of the switching constraint. 
  /// @return Squared norm of the KKT residual in the switching constraint.
  ///
  double KKTError() const;

  ///
  /// @brief Returns l1-norm of the constraint violation, that is, the primal
  /// residual in the switchign constraint. 
  /// @return l1-norm of the constraint violation.
  ///
  double constraintViolation() const;

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
  /// @brief Checks the equivalence of two SwitchingConstraintResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SwitchingConstraintResidual& other) const;

  ///
  /// @brief Checks this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

private:
  Eigen::VectorXd P_full_;
  int dimq_, dimv_, dimi_;

};

} // namespace robotoc 

#include "robotoc/ocp/switching_constraint_residual.hxx"

#endif // ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HPP_ 