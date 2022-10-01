#ifndef ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HPP_
#define ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/sto/sto_cost_function.hpp"
#include "robotoc/sto/sto_constraints.hpp"


namespace robotoc {

///
/// @class SwitchingTimeOptimization
/// @brief The switching time optimization (STO) problem.
///
class SwitchingTimeOptimization {
public:
  ///
  /// @brief Construct the STO problem. 
  /// @param[in] ocp Optimal control problem. 
  ///
  SwitchingTimeOptimization(const OCP& ocp);

  ///
  /// @brief Default Constructor.
  ///
  SwitchingTimeOptimization();

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingTimeOptimization(const SwitchingTimeOptimization&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SwitchingTimeOptimization& operator=(const SwitchingTimeOptimization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingTimeOptimization(SwitchingTimeOptimization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingTimeOptimization& operator=(SwitchingTimeOptimization&&) noexcept = default;

  ///
  /// @brief Reserve the memory for STO problem. 
  /// @param[in] ocp Optimal control problem. 
  ///
  void reserve(const OCP& ocp);

  ///
  /// @brief Sets the regularization paremeter on the STO problem. 
  /// @param[in] sto_reg The regularization paremeter. Must be non-negative.
  ///
  void setRegularization(const double sto_reg);

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in] ocp Optimal control problem.
  ///
  void initConstraints(const OCP& ocp);

  ///
  /// @brief Computes the KKT residual of switching time otpimization problem. 
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTResidual(const OCP& ocp, KKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT system, i.e., the condensed KKT matrix and KKT
  /// residual for Newton's method. 
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTSystem(const OCP& ocp, KKTMatrix& kkt_matrix, 
                        KKTResidual& kkt_residual);

  ///
  /// @return Applies the regularization on the KKT matrix.
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_matrix KKT matrix. 
  ///
  void applyRegularization(const OCP& ocp, KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Compute the Newton direction.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] d Direction. 
  ///
  void computeDirection(const OCP& ocp, const Direction& d);

  // ///
  // /// @brief Returns max primal step size.
  // /// @return max primal step size.
  // /// 
  // double maxPrimalStepSize() const;

  // ///
  // /// @brief Returns max dual step size.
  // /// @return max dual step size.
  // /// 
  // double maxDualStepSize() const;

  ///
  /// @brief Returns the l2-norm of the KKT residual of STO problem.
  ///
  double KKTError() const;

  ///
  /// @brief Returns the l2-norm of the KKT residual of STO problem.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_residual KKT residual. 
  ///
  double KKTError(const OCP& ocp, const KKTResidual& kkt_residual);

  ///
  /// @brief Returns the total value of the cost function.
  ///
  double totalCost() const;

  ///
  /// @brief Integrates the solution (the switching times). The internal 
  /// switching times of contact_sequence are updated. The slack and dual
  /// variables in STO constraints are also updated.
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in, out] d Direction. 
  ///
  void integrateSolution(const OCP& ocp, 
                         std::shared_ptr<ContactSequence>& contact_sequence,
                         const double primal_step_size,
                         const double dual_step_size,
                         const Direction& d) const;




  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in] time_discretization Time discretization. 
  ///
  void initConstraints(const TimeDiscretization& time_discretization);

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in] time_discretization Time discretization. 
  ///
  bool isFeasible(const TimeDiscretization& time_discretization);

  ///
  /// @brief Computes the cost and constraint violations. 
  /// @param[in] time_discretization Time discretization. 
  ///
  void evalSTO(const TimeDiscretization& time_discretization);

  ///
  /// @brief Computes the KKT residual and matrix. 
  /// @param[in] time_discretization Time discretization. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalKKT(const TimeDiscretization& time_discretization, 
               KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  ///
  /// @brief Gets the performance index of the evaluation. 
  /// @return const reference to the performance index.
  ///
  const PerformanceIndex& getEval() const;

  ///
  /// @brief Computes the step sizes via the fraction-to-boundary-rule.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in, out] d Direction. 
  ///
  void computeStepSizes(const TimeDiscretization& time_discretization,
                        Direction& d);

  ///
  /// @brief Gets the maximum primal step size of the fraction-to-boundary-rule.
  /// @return The primal step size of the fraction-to-boundary-rule.
  ///
  double maxPrimalStepSize() const;

  ///
  /// @brief Gets the maximum dual step size of the fraction-to-boundary-rule.
  /// @return The dual step size of the fraction-to-boundary-rule.
  ///
  double maxDualStepSize() const;

  ///
  /// @brief Computes the initial state direction. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in, out] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integrateSolution(const TimeDiscretization& time_discretization, 
                         const double primal_step_size,
                         const double dual_step_size, const Direction& d);

private:
  std::shared_ptr<STOCostFunction> sto_cost_;
  std::shared_ptr<STOConstraints> sto_constraints_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  double sto_reg_, kkt_error_, cost_val_;
  Eigen::VectorXd h_phase_;
  int reserved_num_switches_;
  bool is_sto_enabled_;

  PerformanceIndex performance_index_;
  Eigen::VectorXd lt_, h_, dts_;
  Eigen::MatrixXd Qtt_;
  ConstraintComponentData constraint_data_;
};

} // namespace robotoc

#endif // ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HPP_ 