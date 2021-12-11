#ifndef ROBOTOC_STO_CONSTRAINTS_HPP_
#define ROBOTOC_STO_CONSTRAINTS_HPP_

#include <limits>
#include <cmath>
#include <vector>

#include "Eigen/Core"

#include "robotoc/hybrid/time_discretization.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"
#include "robotoc/ocp/direction.hpp"
#include "robotoc/hybrid/dwell_time_lower_bound.hpp"


namespace robotoc {

///
/// @class STOConstraints
/// @brief Constraintn on the switching times.
///
class STOConstraints {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] max_num_switches Maximum number of switches on the horizon. 
  /// @param[in] min_dt Minimum dwell time. Must be non-negative. The all 
  /// minimum dwell times are set to this value. Default is 
  /// STOConstraints::k_min_dt.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-03.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  STOConstraints(const int max_num_switches, const double min_dt=k_min_dt,
                 const double barrier=1.0e-03, 
                 const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Constructor. 
  /// @param[in] max_num_switches Maximum number of switches on the horizon. 
  /// @param[in] min_dt Minimum dwell times. Each component must be 
  /// non-negative. 
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-03.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  STOConstraints(const int max_num_switches, const std::vector<double>& min_dt,
                 const double barrier=1.0e-03, 
                 const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Default constructor. 
  ///
  STOConstraints();

  ///
  /// @brief Destructor. 
  ///
  ~STOConstraints();

  ///
  /// @brief Default copy constructor. 
  ///
  STOConstraints(const STOConstraints&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  STOConstraints& operator=(const STOConstraints&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  STOConstraints(STOConstraints&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  STOConstraints& operator=(
      STOConstraints&&) noexcept = default;

  ///
  /// @brief Sets the slack variables. 
  /// @param[in] discretization Discretization of the optimal control problem.
  ///
  void setSlack(const TimeDiscretization& discretization);

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] discretization Discretization of the optimal control problem.
  ///
  void evalConstraint(const TimeDiscretization& discretization);

  ///
  /// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
  /// the products of the Jacobian of the constraints and Lagrange multipliers.
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[out] kkt_residual KKT residual.
  ///
  void linearizeConstraints(const TimeDiscretization& discretization,
                            KKTResidual& kkt_residual); 

  ///
  /// @brief Linearizes the constraints (i.e., calls linearizeConstraints())
  /// and condense the slack and dual variables.
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[out] kkt_matrix KKT matrix.
  /// @param[out] kkt_residual KKT residual.
  ///
  void condenseSlackAndDual(const TimeDiscretization& discretization,
                            KKTMatrix& kkt_matrix, 
                            KKTResidual& kkt_residual);

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[in] d Newton direction.
  ///
  void expandSlackAndDual(const TimeDiscretization& discretization, 
                          const Direction& d); 

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the direction of the slack variable.
  /// @return Maximum step size regarding the slack variable.
  ///
  double maxPrimalStepSize() const;

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
  /// @brief Sets the minimum dwell time. 
  /// @param[in] min_dt Minimum dwell time. Must be non-negative. The all 
  /// minimum dwell times are set to this value. Default is 
  /// STOConstraints::k_min_dt.
  ///
  void setMinimumDwellTimes(const double min_dt=k_min_dt);

  ///
  /// @brief Sets the minimum dwell time. 
  /// @param[in] min_dt Minimum dwell times. Each component must be 
  /// non-negative.
  ///
  void setMinimumDwellTimes(const std::vector<double>& min_dt);

  ///
  /// @brief Sets the barrier parameter for all the constraint components.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  ///
  void setBarrier(const double barrier);

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
  double barrier() const;

  ///
  /// @brief Gets the parameter of the fraction-to-boundary-rule. 
  /// @return The parameter of the fraction-to-boundary-rule. 
  ///
  double fractionToBoundaryRule() const;

  ///
  /// @brief Minimum value of the minimum dwell times.  
  ///
  static constexpr double k_min_dt 
      = std::sqrt(std::numeric_limits<double>::epsilon());

private:
  std::vector<DwellTimeLowerBound> dtlb_;
  std::vector<double> min_dt_;
  double k_min_dt_; // this is for Pybind11
  double barrier_, fraction_to_boundary_rule_;
  int max_num_switches_, num_switches_;
  Eigen::VectorXd primal_step_size_, dual_step_size_;

};

} // namespace robotoc

#include "robotoc/hybrid/sto_constraints.hxx"

#endif // ROBOTOC_STO_CONSTRAINTS_HPP_ 