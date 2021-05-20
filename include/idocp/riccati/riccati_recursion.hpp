#ifndef IDOCP_RICCATI_RECURSION_HPP_
#define IDOCP_RICCATI_RECURSION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/utils/aligned_vector.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/split_constrained_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"
#include "idocp/riccati/riccati_factorizer.hpp"
#include "idocp/ocp/ocp.hpp"


namespace idocp {

///
/// @typedef RiccatiFactorization
/// @brief Riccati factorization matices of the LQR subproblem. 
///
using RiccatiFactorization = hybrid_container<SplitRiccatiFactorization, 
                                              SplitRiccatiFactorization, 
                                              SplitConstrainedRiccatiFactorization>;

///
/// @class RiccatiRecursion
/// @brief Riccati recursion solver for hybrid optimal control problems.
///
class RiccatiRecursion {
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
  RiccatiRecursion(const Robot& robot, const int N, const int max_num_impulse, 
                   const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiRecursion();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiRecursion();
 
  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiRecursion(const RiccatiRecursion&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  RiccatiRecursion& operator=(const RiccatiRecursion&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiRecursion(RiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiRecursion& operator=(RiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in, out] factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(const OCP& ocp, KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual, 
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
  static void computeInitialStateDirection(const aligned_vector<Robot>& robots, 
                                           const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v, 
                                           const KKTMatrix& kkt_matrix, 
                                           const Solution& s, Direction& d);

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] d Direction. d[0].dx must be computed using 
  /// computeInitialStateDirection.
  ///
  void forwardRiccatiRecursion(const OCP& ocp, const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual, 
                               Direction& d) const;

  ///
  /// @brief Compute the Newton direction in parallel from the Riccati 
  /// factorization factorized by 
  /// RiccatiRecursion::backwardRiccatiRecursion() and 
  /// RiccatiRecursion::forwardRiccatiRecursion().
  /// @param[in] ocp Optimal control problem.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] factorization Riccati factorization. 
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  void computeDirection(OCP& ocp, aligned_vector<Robot>& robots, 
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
  int nthreads_, N_, N_all_;
  RiccatiFactorizer factorizer_;
  hybrid_container<LQRPolicy> lqr_policy_;
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;

};

} // namespace idocp

#endif // IDOCP_RICCATI_RECURSION_HPP_ 