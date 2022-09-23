#ifndef ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HPP_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HPP_

#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"


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
  /// @brief Default destructor. 
  ///
  ~SwitchingConstraintResidual() = default;

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
  /// @brief Sets the dimension of the switching constraint.
  /// @param[in] dims The dimension of the switching constraint. Must be non-negative.
  ///
  void setDimension(const int dims);

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
  /// @brief Returns the lp norm of the constraint violation, that is,
  /// the primal residual in the switching constraint. Default norm is l1-norm.
  /// You can specify l-infty norm by passing Eigen::Infinity as the 
  /// template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the constraint violation.
  ///
  template <int p=1>
  double constraintViolation() const;

  template <int p=1>
  double primalFeasibility() const {
    if (P().size() > 0)
      P().template lpNorm<p>();
    else
      return 0.0;
  }

  template <int p=1>
  static double dualFeasibility() {
    return 0.0;
  }

  ///
  /// @brief Sets the split KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the switching constraint.
  /// @return Dimension of the switching constraint.
  ///
  int dims() const;

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

  ///
  /// @brief Displays the switching constraint residual onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const SwitchingConstraintResidual& sc_residual);

private:
  Eigen::VectorXd P_full_;
  int dimq_, dimv_, dims_;

};

} // namespace robotoc 

#include "robotoc/core/switching_constraint_residual.hxx"

#endif // ROBOTOC_SWITCHING_CONSTRAINT_RESIDUAL_HPP_ 