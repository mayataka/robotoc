#ifndef IDOCP_RICCATI_SOLVER_HPP_
#define IDOCP_RICCATI_SOLVER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/ocp/riccati_direction_calculator.hpp"


namespace idocp {

///
/// @class RiccatiSolver
/// @brief Riccati solver of the optimal control problem. 
///
class RiccatiSolver {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  RiccatiSolver(const Robot& robot, const int N, const int max_num_impulse, 
                const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiSolver();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiSolver();

  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiSolver(const RiccatiSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  RiccatiSolver& operator=(const RiccatiSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiSolver(RiccatiSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiSolver& operator=(RiccatiSolver&&) noexcept = default;

  ///
  /// @brief Computes the Newton direction using Riccati recursion. Before 
  /// calling this function, call OCPLinearizer::linearizeOCP() to compute
  /// kkt_matrix and kkt_residual.
  /// @param[in] ocp OCP. 
  /// @param[in] robots Robot models. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration vector. 
  /// @param[in] v Initial generalized velocity vector. 
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in] sampling_period Sampling period. Default is 0.
  /// 
  template <bool UseContinuationMethod=false>
  void computeNewtonDirection(OCP& ocp, std::vector<Robot>& robots, 
                              const ContactSequence& contact_sequence,
                              const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, const Solution& s, 
                              Direction& d, KKTMatrix& kkt_matrix, 
                              KKTResidual& kkt_residual,
                              const double sampling_period=0);

  ///
  /// @brief Computes the maximum step size computed from the primal variables. 
  /// Before call this function, call RiccatiSolver::computeNewtonDirection().
  /// @return Maximum primal step size
  /// 
  double maxPrimalStepSize() const;

  ///
  /// @brief Computes the maximum step size computed from the dual variables. 
  /// Before call this function, call RiccatiSolver::computeNewtonDirection().
  /// @return maximum dual step size
  /// 
  double maxDualStepSize() const;

  ///
  /// @brief Getter of the state feedback gain of the LQR subproblem at the 
  /// specified time stage. 
  /// @param[in] Kq The state feedback gain with respect to the configuration. 
  /// @param[in] Kv The state feedback gain with respect to the velocity. 
  /// @param[in] time_stage Time stage of interested. 
  ///
  void getStateFeedbackGain(const int time_stage, Eigen::MatrixXd& Kq, 
                            Eigen::MatrixXd& Kv) const;

private:
  RiccatiRecursion riccati_recursion_;
  RiccatiFactorizer riccati_factorizer_;
  RiccatiFactorization riccati_factorization_;
  StateConstraintRiccatiFactorizer constraint_factorizer_;
  StateConstraintRiccatiFactorization constraint_factorization_;
  RiccatiDirectionCalculator direction_calculator_;

};

} // namespace idocp 

#include "idocp/ocp/riccati_solver.hxx"

#endif // IDOCP_RICCATI_SOLVER_HPP_ 