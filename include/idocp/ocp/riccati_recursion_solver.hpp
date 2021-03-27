#ifndef IDOCP_RICCATI_RECURSION_SOLVER_HPP_
#define IDOCP_RICCATI_RECURSION_SOLVER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"

namespace idocp {

///
/// @class RiccatiRecursionSolver
/// @brief Riccati recursion solver for hybrid optimal control problems.
///
class RiccatiRecursionSolver {
public:
  ///
  /// @brief Construct a Riccati recursion solver.
  /// @param[in] robot Robot model. 
  /// @param[in] N Number of discretization of the horizon. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  RiccatiRecursionSolver(const Robot& robot, const int N, 
                         const int max_num_impulse, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiRecursionSolver();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiRecursionSolver();
 
  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiRecursionSolver(const RiccatiRecursionSolver&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  RiccatiRecursionSolver& operator=(const RiccatiRecursionSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiRecursionSolver(RiccatiRecursionSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiRecursionSolver& operator=(RiccatiRecursionSolver&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in] jac Jacobian of the switching constraints. 
  /// @param[in, out] factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(const OCP& ocp, KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual, 
                                const StateConstraintJacobian& jac,
                                RiccatiFactorization& factorization);

  ///
  /// @brief Computes initial state direction.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  static void computeInitialStateDirection(const std::vector<Robot>& robots, 
                                           const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v, 
                                           const KKTMatrix& kkt_matrix, 
                                           const Solution& s, Direction& d);

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] d Direction. d[0].dx() must be computed using 
  /// computeInitialStateDirection.
  ///
  void forwardRiccatiRecursion(const OCP& ocp, const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual, 
                               Direction& d) const;

  ///
  /// @brief Compute the Newton direction in parallel from the Riccati 
  /// factorization factorized by 
  /// RiccatiRecursionSolver::backwardRiccatiRecursion() and 
  /// RiccatiRecursionSolver::forwardRiccatiRecursion().
  /// @param[in] ocp Optimal control problem.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] factorization Riccati factorization. 
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  void computeDirection(OCP& ocp, std::vector<Robot>& robots, 
                        const RiccatiFactorization& factorization, 
                        const Solution& s, Direction& d);

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
  /// @brief Gets of the state feedback gain of the LQR subproblem of the 
  /// specified time stage. 
  /// @param[in] time_stage Time stage of interested. 
  /// @param[in, out] Kq The state feedback gain with respect to the configuration. 
  /// @param[in, out] Kv The state feedback gain with respect to the velocity. 
  ///
  void getStateFeedbackGain(const int time_stage, Eigen::MatrixXd& Kq, 
                            Eigen::MatrixXd& Kv) const;

private:
  int nthreads_, N_all_;
  RiccatiFactorizer factorizer_;
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;

};

} // namespace idocp

#endif // IDOCP_RICCATI_RECURSION_SOLVER_HPP_ 