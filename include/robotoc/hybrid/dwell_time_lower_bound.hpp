#ifndef ROBOTOC_DWELL_TIME_LOWER_BOUND_HPP_
#define ROBOTOC_DWELL_TIME_LOWER_BOUND_HPP_

#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class DwellTimeLowerBound
/// @brief Constraint on the lower limit on the dwell time of the each interval.
///
class DwellTimeLowerBound {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. 
  ///
  DwellTimeLowerBound(const double barrier_param, 
                      const double fraction_to_boundary_rule);

  ///
  /// @brief Default constructor. 
  ///
  DwellTimeLowerBound();

  ///
  /// @brief Destructor. 
  ///
  ~DwellTimeLowerBound();

  ///
  /// @brief Default copy constructor. 
  ///
  DwellTimeLowerBound(const DwellTimeLowerBound&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  DwellTimeLowerBound& operator=(const DwellTimeLowerBound&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  DwellTimeLowerBound(DwellTimeLowerBound&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  DwellTimeLowerBound& operator=(DwellTimeLowerBound&&) noexcept = default;

  ///
  /// @brief Sets the slack variables. 
  /// @param[in] min_dt Lower bound of the dwell time.
  /// @param[in] ts1 The former switching time of the horizon.
  /// @param[in] ts2 The latter switching time of the horizon.
  ///
  void setSlack(const double min_dt, const double ts1, const double ts2);

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] min_dt Lower bound of the dwell time.
  /// @param[in] ts1 The former switching time of the horizon.
  /// @param[in] ts2 The latter switching time of the horizon.
  ///
  void evalConstraint(const double min_dt, const double ts1, const double ts2);

  ///
  /// @brief Computes the derivatives of the priaml residual, i.e., the 
  /// Jacobian of the inequality constraint, and add the product of the 
  /// Jacobian and the dual variable to the KKT residual. This function is 
  /// always called just after evalConstraint().
  /// @param[out] kkt_residual1 Split KKT residual of the former switch.
  /// @param[out] kkt_residual2 Split KKT residual of the latter switch.
  ///
  void evalDerivatives_lub(SplitKKTResidual& kkt_residual1, 
                           SplitKKTResidual& kkt_residual2) const; 

  ///
  /// @brief Computes the derivatives of the priaml residual, i.e., the 
  /// Jacobian of the inequality constraint, and add the product of the 
  /// Jacobian and the dual variable to the KKT residual. This function is 
  /// always called just after evalConstraint().
  /// @param[out] kkt_residual1 Split KKT residual of the former switch.
  ///
  void evalDerivatives_ub(SplitKKTResidual& kkt_residual1) const; 

  ///
  /// @brief Computes the derivatives of the priaml residual, i.e., the 
  /// Jacobian of the inequality constraint, and add the product of the 
  /// Jacobian and the dual variable to the KKT residual. This function is 
  /// always called just after evalConstraint().
  /// @param[out] kkt_residual2 Split KKT residual of the latter switch.
  ///
  void evalDerivatives_lb(SplitKKTResidual& kkt_residual2) const; 

  ///
  /// @brief Condenses the slack and dual variables, i.e., factorizes the  
  /// condensed Hessians and KKT residuals. This function is always called 
  /// just after evalDerivatives().
  /// @param[out] kkt_matrix1 Split KKT matrix of the former switch.
  /// @param[out] kkt_residual1 Split KKT residual of the former switch.
  /// @param[out] kkt_matrix2 Split KKT matrix of the latter switch.
  /// @param[out] kkt_residual2 Split KKT residual of the latter switch.
  ///
  void condenseSlackAndDual_lub(SplitKKTMatrix& kkt_matrix1, 
                                SplitKKTResidual& kkt_residual1,
                                SplitKKTMatrix& kkt_matrix2, 
                                SplitKKTResidual& kkt_residual2) const;

  ///
  /// @brief Condenses the slack and dual variables, i.e., factorizes the  
  /// condensed Hessians and KKT residuals. This function is always called 
  /// just after evalDerivatives().
  /// @param[out] kkt_matrix1 Split KKT matrix of the former switch.
  /// @param[out] kkt_residual1 Split KKT residual of the latter switch.
  ///
  void condenseSlackAndDual_ub(SplitKKTMatrix& kkt_matrix1, 
                               SplitKKTResidual& kkt_residual1) const;

  ///
  /// @brief Condenses the slack and dual variables, i.e., factorizes the  
  /// condensed Hessians and KKT residuals. This function is always called 
  /// just after evalDerivatives().
  /// @param[out] kkt_matrix2 Split KKT matrix of the former switch.
  /// @param[out] kkt_residual2 Split KKT residual of the latter switch.
  ///
  void condenseSlackAndDual_lb(SplitKKTMatrix& kkt_matrix2, 
                               SplitKKTResidual& kkt_residual2) const;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] dts1 Direction of the former switching time.
  /// @param[in] dts2 Direction of the latter switching time.
  ///
  void expandSlackAndDual_lub(const double dts1, const double dts2); 

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] dts1 Direction of the former switching time.
  ///
  void expandSlackAndDual_ub(const double dts1); 

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] dts2 Direction of the latter switching time.
  ///
  void expandSlackAndDual_lb(const double dts2); 

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the direction of the slack variable.
  /// @return Maximum step size regarding the slack variable.
  ///
  double maxSlackStepSize() const;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the direction of the dual variable.
  /// @param[in] data Constraint data. 
  /// @return Maximum step size regarding the dual variable.
  ///
  double maxDualStepSize() const;

  ///
  /// @brief Updates the slack variable according to the step size.
  /// @param[in] step_size Step size. 
  ///
  void updateSlack(const double step_size);

  ///
  /// @brief Updates the dual variable according to the step size.
  /// @param[in] step_size Step size. 
  ///
  void updateDual(const double step_size);

  ///
  /// @brief Returns the sum of the squared norm of the KKT error 
  /// (primal residual and complementary slackness) of all the constraints. 
  /// @return KKT error. 
  ///
  double KKTError() const;

  ///
  /// @brief Sets the barrier parameter for all the constraint components.
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  ///
  void setBarrierParam(const double barrier_param);

  ///
  /// @brief Sets the parameter of the fraction-to-boundary-rule for all the 
  /// constraint components.
  /// @param[in] fraction_to_boundary_rule Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995.
  ///
  void setFractionToBoundaryRule(const double fraction_to_boundary_rule);

  ///
  /// @brief Gets the barrier parameter.
  /// @return Barrier parameter. 
  ///
  double getBarrierParam() const;

  ///
  /// @brief Gets the parameter of the fraction-to-boundary-rule. 
  /// @return The parameter of the fraction-to-boundary-rule. 
  ///
  double getFractionToBoundaryRule() const;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const DwellTimeLowerBound& dtlb);

private:
  double barrier_, fraction_to_boundary_rule_, 
         slack_, dual_, residual_, cmpl_, dslack_, ddual_, log_barrier_;

};

} // namespace robotoc

#include "robotoc/hybrid/dwell_time_lower_bound.hxx"

#endif // ROBOTOC_DWELL_TIME_LOWER_BOUND_HPP_ 