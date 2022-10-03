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
/// @brief Constraintn on the switching times.
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

  ConstraintComponentData createConstraintsData(
      const TimeDiscretization& time_discretization) const;

  bool isFeasible(const TimeDiscretization& time_discretization,
                  ConstraintComponentData& data) const;

  void setSlackAndDual(const TimeDiscretization& time_discretization,
                       ConstraintComponentData& data) const;

  void evalConstraint(const TimeDiscretization& time_discretization, 
                      ConstraintComponentData& data) const;

  void linearizeConstraints(const TimeDiscretization& time_discretization, 
                            ConstraintComponentData& data, 
                            Eigen::VectorXd& lt) const;

  void condenseSlackAndDual(ConstraintComponentData& data, Eigen::VectorXd& lt,
                            Eigen::MatrixXd& Qtt) const;

  void expandSlackAndDual(ConstraintComponentData& data, Eigen::VectorXd& dts) const;

  double maxSlackStepSize(const ConstraintComponentData& data) const;

  double maxDualStepSize(const ConstraintComponentData& data) const;

  ///
  /// @brief Updates the slack variable according to the step size.
  /// @param[in] step_size Step size. 
  ///
  void updateSlack(ConstraintComponentData& data, const double step_size) const;

  ///
  /// @brief Updates the dual variable according to the step size.
  /// @param[in] step_size Step size. 
  ///
  void updateDual(ConstraintComponentData& data, const double step_size) const;

  static void computeDwellTimes(const TimeDiscretization& time_discretization,
                                Eigen::VectorXd& dwell_times);

private:
  double barrier_, fraction_to_boundary_rule_, eps_;
  int num_switches_;

  Eigen::VectorXd primal_step_size_, dual_step_size_;
  Eigen::VectorXd minimum_dwell_times_;

};

} // namespace robotoc

#endif // ROBOTOC_STO_CONSTRAINTS_HPP_ 