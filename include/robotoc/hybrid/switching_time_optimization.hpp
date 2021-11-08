#ifndef ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HPP_
#define ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/hybrid/hybrid_ocp_discretization.hpp"
#include "robotoc/hybrid/sto_cost_function.hpp"
#include "robotoc/hybrid/sto_constraints.hpp"
#include "robotoc/hybrid/sto_regularization.hpp"


namespace robotoc {

///
/// @class SwitchingTimeOptimization
/// @brief The switching time optimization (STO) problem.
///
class SwitchingTimeOptimization {
public:
  ///
  /// @brief Construct the STO problem. 
  /// @param[in] sto_cost Shared ptr to the STO cost function.
  /// @param[in] sto_constraints Shared ptr to the STO constraints.
  /// @param[in] max_num_impulse_events Maximum number of the impulse on the 
  /// horizon. Must be non-negative. 
  ///
  SwitchingTimeOptimization(
      const std::shared_ptr<STOCostFunction>& sto_cost, 
      const std::shared_ptr<STOConstraints>& sto_constraints, 
      const int max_num_impulse_events);

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
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in] ocp Optimal control problem.
  ///
  void initConstraints(const OCP& ocp) const;

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
  /// @param[in] kkt_error KKT error. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  ///
  void applyRegularization(const OCP& ocp, const double kkt_error, 
                           KKTMatrix& kkt_matrix) const;

  ///
  /// @brief Compute the Newton direction.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] d Direction. 
  ///
  void computeDirection(const OCP& ocp, const Direction& d);

  ///
  /// @brief Returns max primal step size.
  /// @return max primal step size.
  /// 
  double maxPrimalStepSize() const;

  ///
  /// @brief Returns max dual step size.
  /// @return max dual step size.
  /// 
  double maxDualStepSize() const;

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
  /// @brief Set the regularization for the STO problem
  /// @param[in] sto_reg Regularization for the STO problem.
  ///
  void setSTORegularization(const STORegularization& sto_reg);

  ///
  /// @brief Sets regularization method.
  /// @param[in] reg_type regularization type. 
  /// @param[in] w weight parameter (a.k.a. scaling parameter) of the 
  /// regularization. Must be non-negative.
  ///
  void setSTORegularization(const STORegularizationType& reg_type,
                            const double w);

private:
  std::shared_ptr<STOCostFunction> sto_cost_;
  std::shared_ptr<STOConstraints> sto_constraints_;
  STORegularization sto_reg_;
  int max_num_impulse_events_;
  double kkt_error_, cost_val_;
  Eigen::VectorXd h_phase_;
  bool is_sto_enabled_;
};

} // namespace robotoc

#include "robotoc/hybrid/switching_time_optimization.hxx"

#endif // ROBOTOC_SWITCHING_TIME_OPTIMIZATION_HPP_ 