#ifndef ROBOTOC_UNCONSTR_DIRECT_MULTIPLE_SHOOTING_HPP_
#define ROBOTOC_UNCONSTR_DIRECT_MULTIPLE_SHOOTING_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/core/solution.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/grid_info.hpp"
#include "robotoc/unconstr/unconstr_ocp_data.hpp"
#include "robotoc/unconstr/unconstr_intermediate_stage.hpp"
#include "robotoc/unconstr/unconstr_terminal_stage.hpp"


namespace robotoc {

///
/// @class UnconstrDirectMultipleShooting
/// @brief Direct multiple shooting method of hybrid optimal control problems
/// of unconstrained rigid-body systems. 
///
class UnconstrDirectMultipleShooting {
public:
  ///
  /// @brief Construct the direct multiple shooting method.
  /// @param[in] ocp Optimal control problem. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  UnconstrDirectMultipleShooting(const OCP& ocp, const int nthreads=1);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrDirectMultipleShooting();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrDirectMultipleShooting() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrDirectMultipleShooting(
      const UnconstrDirectMultipleShooting&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnconstrDirectMultipleShooting& operator=(
      const UnconstrDirectMultipleShooting&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrDirectMultipleShooting(
      UnconstrDirectMultipleShooting&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrDirectMultipleShooting& operator=(
      UnconstrDirectMultipleShooting&&) noexcept = default;

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] s Solution. 
  ///
  void initConstraints(aligned_vector<Robot>& robots,
                       const std::vector<GridInfo>& time_discretization, 
                       const Solution& s);

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] s Solution. 
  ///
  bool isFeasible(aligned_vector<Robot>& robots, 
                  const std::vector<GridInfo>& time_discretization, 
                  const Solution& s);

  ///
  /// @brief Computes the cost and constraint violations. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalOCP(aligned_vector<Robot>& robots, 
               const std::vector<GridInfo>& time_discretization, 
               const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
               const Solution& s, KKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual and matrix. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalKKT(aligned_vector<Robot>& robots, 
               const std::vector<GridInfo>& time_discretization, 
               const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
               const Solution& s, KKTMatrix& kkt_matrix, 
               KKTResidual& kkt_residual);

  ///
  /// @brief Computes the initial state direction. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  static void computeInitialStateDirection(const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v, 
                                           const Solution& s, Direction& d);

  ///
  /// @brief Gets the performance index of the evaluation. 
  /// @return const reference to the performance index.
  ///
  const PerformanceIndex& getEval() const;

  ///
  /// @brief Computes the step sizes via the fraction-to-boundary-rule.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in, out] d Direction. 
  ///
  void computeStepSizes(const std::vector<GridInfo>& time_discretization,
                        KKTMatrix& kkt_matrix, KKTResidual& kkt_residual,
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
  /// @brief Integrates the solution. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in, out] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integrateSolution(const aligned_vector<Robot>& robots,
                         const std::vector<GridInfo>& time_discretization, 
                         const double primal_step_size,
                         const double dual_step_size,
                         Direction& d, Solution& s);

private:
  int nthreads_;
  UnconstrIntermediateStage intermediate_stage_;
  UnconstrTerminalStage terminal_stage_;
  aligned_vector<UnconstrOCPData> data_;
  PerformanceIndex performance_index_; 
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;
};

} // namespace robotoc 

#endif // ROBOTOC_UNCONSTR_DIRECT_MULTIPLE_SHOOTING_HPP_