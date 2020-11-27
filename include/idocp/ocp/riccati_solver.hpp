#ifndef IDOCP_RICCATI_SOLVER_HPP_
#define IDOCP_RICCATI_SOLVER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/ocp/ocp_direction_calculator.hpp"


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
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  RiccatiSolver(const Robot& robot, const double T, const int N, 
                const int max_num_impulse, const int num_proc);

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
  /// @brief Sets the dimensions of the pure-state equality constraints from
  /// contact sequence. 
  /// @param[in] contact_sequence Contact sequence. 
  /// 
  void setConstraintDimensions(const ContactSequence& contact_sequence);

  ///
  /// @brief Computes the Newton direction using Riccati recursion. Call after
  /// calling OCPLinearizer::linearizeOCP() before calling this function.
  /// @param[in] split_ocps Split OCPs. 
  /// @param[in] robots Robot models. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration vector. 
  /// @param[in] v Initial generalized velocity vector. 
  /// @param[in] s Splits solutions. 
  /// @param[in] d Splits directions. 
  /// @param[in] kkt_matrix KKT matrices. 
  /// @param[in] kkt_residual KKT residuals. 
  /// 
  void computeNewtonDirection(HybridOCP& split_ocps, std::vector<Robot>& robots, 
                              const ContactSequence& contact_sequence,
                              const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, const HybridSolution& s, 
                              HybridDirection& d, HybridKKTMatrix& kkt_matrix, 
                              HybridKKTResidual& kkt_residual);

  ///
  /// @brief Computes the maximum step size computed from the primal variables. 
  /// Before call this function, call RiccatiSolver::computeNewtonDirection().
  /// @param[in] contact_sequence Contact sequence. 
  /// @return Maximum primal step size
  /// 
  double maxPrimalStepSize(const ContactSequence& contact_sequence) const;

  ///
  /// @brief Computes the maximum step size computed from the dual variables. 
  /// Before call this function, call RiccatiSolver::computeNewtonDirection().
  /// @param[in] contact_sequence Contact sequence. 
  /// @return maximum dual step size
  /// 
  double maxDualStepSize(const ContactSequence& contact_sequence) const;

private:
  RiccatiRecursion riccati_recursion_;
  HybridRiccatiFactorizer riccati_factorizer_;
  HybridRiccatiFactorization riccati_factorization_;
  StateConstraintRiccatiFactorizer constraint_factorizer_;
  StateConstraintRiccatiFactorization constraint_factorization_;
  OCPDirectionCalculator ocp_direction_calculator_;

};

} // namespace idocp 

#endif // IDOCP_RICCATI_SOLVER_HPP_ 