#ifndef IDOCP_SWITCHING_TIME_CONSTRAINTS_HPP_
#define IDOCP_SWITCHING_TIME_CONSTRAINTS_HPP_

#include <limits>
#include <cmath>
#include <vector>

#include "Eigen/Core"

#include "idocp/hybrid/hybrid_time_discretization.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/direction.hpp"
#include "idocp/hybrid/dwell_time_lower_bound.hpp"


namespace idocp {

///
/// @class SwitchingTimeConstraints
/// @brief Constraintn on the switching times.
///
class SwitchingTimeConstraints {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] max_num_switches Maximum number of switches on the horizon. 
  /// @param[in] min_dt Minimum dwell-time. Default is 
  /// std::sqrt(std::numeric_limits<double>::epsilon()).
  /// @param[in] min_dt0 Minimum time margin between the initial time of the 
  /// horizon and initial switching time. Default is 
  /// std::sqrt(std::numeric_limits<double>::epsilon()).
  /// @param[in] min_dtf Minimum time margin between the terminal time of the 
  /// horizon and last switching time. Default is 
  /// std::sqrt(std::numeric_limits<double>::epsilon()).
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-04.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  SwitchingTimeConstraints(
      const int max_num_switches,
      const double min_dt=std::sqrt(std::numeric_limits<double>::epsilon()),
      const double min_dt0=std::sqrt(std::numeric_limits<double>::epsilon()),
      const double min_dtf=std::sqrt(std::numeric_limits<double>::epsilon()),
      const double barrier=1.0e-04, 
      const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Default constructor. 
  ///
  SwitchingTimeConstraints();

  ///
  /// @brief Destructor. 
  ///
  ~SwitchingTimeConstraints();

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingTimeConstraints(const SwitchingTimeConstraints&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SwitchingTimeConstraints& operator=(const SwitchingTimeConstraints&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingTimeConstraints(SwitchingTimeConstraints&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingTimeConstraints& operator=(
      SwitchingTimeConstraints&&) noexcept = default;

  ///
  /// @brief Sets the slack variables. 
  /// @param[in] discretization Discretization of the optimal control problem.
  ///
  void setSlack(const HybridTimeDiscretization& discretization);

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] discretization Discretization of the optimal control problem.
  ///
  void evalConstraint(const HybridTimeDiscretization& discretization);

  ///
  /// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
  /// the products of the Jacobian of the constraints and Lagrange multipliers.
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[out] kkt_residual KKT residual.
  ///
  void linearizeConstraints(const HybridTimeDiscretization& discretization,
                            KKTResidual& kkt_residual); 

  ///
  /// @brief Linearizes the constraints (i.e., calls linearizeConstraints())
  /// and condense the slack and dual variables.
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[out] kkt_matrix KKT matrix.
  /// @param[out] kkt_residual KKT residual.
  ///
  void condenseSlackAndDual(const HybridTimeDiscretization& discretization,
                            KKTMatrix& kkt_matrix, 
                            KKTResidual& kkt_residual);

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] discretization Discretization of the optimal control problem.
  /// @param[in] d Newton direction.
  ///
  void expandSlackAndDual(const HybridTimeDiscretization& discretization, 
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

private:
  std::vector<DwellTimeLowerBound> dtlb_;
  double min_dt_, min_dt0_, min_dtf_;
  int num_switches_;
  Eigen::VectorXd primal_step_size_, dual_step_size_;

};

} // namespace idocp

#include "idocp/hybrid/switching_time_constraints.hxx"

#endif // IDOCP_SWITCHING_TIME_CONSTRAINTS_HPP_ 