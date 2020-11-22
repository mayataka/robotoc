#ifndef IDOCP_OCP_RICCATI_SOLVER_HPP_
#define IDOCP_OCP_RICCATI_SOLVER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/ocp/riccati_recursion.hpp"
#include "idocp/ocp/ocp_direction_calculator.hpp"


namespace idocp {

///
/// @class OCPRiccatiSolver
/// @brief Linearize of the optimal control problem. 
///
class OCPRiccatiSolver {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPRiccatiSolver(const Robot& robot, const double T, const int N, 
                   const int max_num_impulse, const int num_proc);

  ///
  /// @brief Default constructor. 
  ///
  OCPRiccatiSolver();

  ///
  /// @brief Destructor. 
  ///
  ~OCPRiccatiSolver();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPRiccatiSolver(const OCPRiccatiSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPRiccatiSolver& operator=(const OCPRiccatiSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPRiccatiSolver(OCPRiccatiSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPRiccatiSolver& operator=(OCPRiccatiSolver&&) noexcept = default;

  void computeDirection(
      HybridOCP& split_ocps, std::vector<Robot>& robots, 
      const ContactSequence& contact_sequence,
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
      const HybridSolution& s, HybridDirection& d,
      HybridKKTMatrix& kkt_matrix, HybridKKTResidual& kkt_residual);

  double maxPrimalStepSize(const ContactSequence& contact_sequence) const;

  double maxDualStepSize(const ContactSequence& contact_sequence) const;

private:
  RiccatiRecursion riccati_recursion_;
  HybridRiccatiFactorization riccati_factorization_;
  StateConstraintRiccatiFactorizer constraint_factorizer_;
  std::vector<StateConstraintRiccatiFactorization> constraint_factorization_;
  OCPDirectionCalculator ocp_direction_calculator_;

};

} // namespace idocp 

#include "idocp/ocp/ocp_riccati_solver.hxx"

#endif // IDOCP_OCP_RICCATI_SOLVER_HPP_ 