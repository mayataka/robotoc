#ifndef ROBOTOC_STO_CONSTRAINTS_HPP_
#define ROBOTOC_STO_CONSTRAINTS_HPP_

#include <limits>
#include <cmath>
#include <vector>

#include "Eigen/Core"

#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"


namespace robotoc {

///
/// @class STOConstraints
/// @brief Constraints of the switching time optimization problems.
///
class STOConstraints {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] num_switches The number of switches.
  /// @param[in] minimum_dwell_time Minimum dwell time. Must be non-negative. 
  /// The all minimum dwell times are set to this value. Default is 
  /// std::sqrt(std::numeric_limits<double>::epsilon()).
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-03.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  STOConstraints(const int num_switches, 
                 const double minimum_dwell_time=std::sqrt(std::numeric_limits<double>::epsilon()),
                 const double barrier_param=1.0e-03, 
                 const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Constructor. 
  /// @param[in] minimum_dwell_times Minimum dwell times. Each component must be 
  /// non-negative. 
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-03.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  STOConstraints(const std::vector<double>& minimum_dwell_times,
                 const double barrier_param=1.0e-03, 
                 const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Constructor. 
  /// @param[in] minimum_dwell_times Minimum dwell times. Each component must be 
  /// non-negative. 
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-03.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  STOConstraints(const Eigen::VectorXd& minimum_dwell_times,
                 const double barrier_param=1.0e-03, 
                 const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Default constructor. 
  ///
  STOConstraints();

  ///
  /// @brief Default destructor. 
  ///
  ~STOConstraints() = default;

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
  STOConstraints& operator=(STOConstraints&&) noexcept = default;

  ///
  /// @brief Sets the minimum dwell times. 
  /// @param[in] minimum_dwell_time Minimum dwell time. Must be non-negative. 
  /// The all minimum dwell times are set to this value. Default is 
  /// std::sqrt(std::numeric_limits<double>::epsilon()).
  ///
  void setMinimumDwellTimes(
      const double minimum_dwell_time=std::sqrt(std::numeric_limits<double>::epsilon()));

  ///
  /// @brief Sets the minimum dwell times. 
  /// @param[in] minimum_dwell_times Minimum dwell times. Each component must be 
  /// non-negative. 
  ///
  void setMinimumDwellTimes(const std::vector<double>& minimum_dwell_times);

  ///
  /// @brief Sets the minimum dwell times. 
  /// @param[in] minimum_dwell_times Minimum dwell times. Each component must be 
  /// non-negative. 
  ///
  void setMinimumDwellTimes(const Eigen::VectorXd& minimum_dwell_times);

  ///
  /// @brief Gets the minimum dwell times. 
  /// @return const reference to the minimum dwell times. 
  ///
  const Eigen::VectorXd& getMinimumDwellTimes() const;

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

  ///
  /// @brief Creates the data for this constraints. 
  /// @param[in] time_discretization Time discretization.
  /// @return The constraint data for this constraints.
  ///
  ConstraintComponentData createConstraintsData(
      const TimeDiscretization& time_discretization) const;

  ///
  /// @brief Checks whether the current split solution s is feasible or not. 
  /// @param[in] time_discretization Time discretization.
  /// @param[in, out] data Constraints data. 
  /// @return true if s is feasible. false if not.
  ///
  bool isFeasible(const TimeDiscretization& time_discretization,
                  ConstraintComponentData& data) const;

  ///
  /// @brief Sets the slack and dual variables of constraints. 
  /// @param[in] time_discretization Time discretization.
  /// @param[in, out] data Constraints data. 
  ///
  void setSlackAndDual(const TimeDiscretization& time_discretization,
                       ConstraintComponentData& data) const;

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] time_discretization Time discretization.
  /// @param[in, out] data Constraints data. 
  ///
  void evalConstraint(const TimeDiscretization& time_discretization, 
                      ConstraintComponentData& data) const;

  ///
  /// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
  /// the products of the Jacobian of the constraints and Lagrange multipliers.
  /// @param[in] time_discretization Time discretization.
  /// @param[in, out] data Constraints data. 
  /// @param[in, out] lt The derivatives of the Lagrangian w.r.t. the switching
  /// times. 
  ///
  void linearizeConstraints(const TimeDiscretization& time_discretization, 
                            ConstraintComponentData& data, 
                            Eigen::VectorXd& lt) const;

  ///
  /// @brief Condenses the slack and dual variables. linearizeConstraints() must 
  /// be called before this function.
  /// @param[in, out] data Constraints data. 
  /// @param[in, out] lt The derivatives of the Lagrangian w.r.t. the switching
  /// times. 
  /// @param[in, out] Qtt The Hessian of the Lagrangian w.r.t. the switching
  /// times. 
  ///
  void condenseSlackAndDual(ConstraintComponentData& data, Eigen::VectorXd& lt,
                            Eigen::MatrixXd& Qtt) const;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in, out] data Constraints data. 
  /// @param[in, out] dts The direction of the switching times.
  ///
  void expandSlackAndDual(ConstraintComponentData& data, Eigen::VectorXd& dts) const;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the directions of the slack variables.
  /// @param[in] data Constraints data.
  /// @return Maximum step size regarding the slack variables.
  ///
  double maxSlackStepSize(const ConstraintComponentData& data) const;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the directions of the dual variables.
  /// @param[in] data Constraints data.
  /// @return Maximum step size regarding the dual variables.
  ///
  double maxDualStepSize(const ConstraintComponentData& data) const;

  ///
  /// @brief Updates the slack variables according to step_size.
  /// @param[in, out] data Constraints data. 
  /// @param[in] step_size Step size. 
  ///
  void updateSlack(ConstraintComponentData& data, const double step_size) const;

  ///
  /// @brief Updates the dual variables according to step_size.
  /// @param[in, out] data Constraints data. 
  /// @param[in] step_size Step size. 
  ///
  void updateDual(ConstraintComponentData& data, const double step_size) const;

  ///
  /// @brief Computes the dwell times from the given time discretization.
  /// @param[in] time_discretization Time discretization.
  /// @param[in, out] dwell_times Dwell times. 
  ///
  static void computeDwellTimes(const TimeDiscretization& time_discretization,
                                Eigen::VectorXd& dwell_times);

private:
  double barrier_, fraction_to_boundary_rule_;
  int num_switches_;

  Eigen::VectorXd primal_step_size_, dual_step_size_;
  Eigen::VectorXd minimum_dwell_times_;

};

} // namespace robotoc

#endif // ROBOTOC_STO_CONSTRAINTS_HPP_ 