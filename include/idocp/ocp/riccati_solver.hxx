#ifndef IDOCP_RICCATI_SOLVER_HXX_
#define IDOCP_RICCATI_SOLVER_HXX_

#include "idocp/ocp/riccati_solver.hpp"

namespace idocp {

template <bool UseContinuationMethod>
inline void RiccatiSolver::computeNewtonDirection(
    OCP& ocp, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, Direction& d, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    const double sampling_period) {
  riccati_recursion_.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, 
                                                      riccati_factorization_);
  riccati_recursion_.backwardRiccatiRecursion(riccati_factorizer_, 
                                              ocp.discrete(), kkt_matrix, 
                                              kkt_residual, 
                                              riccati_factorization_);
  constraint_factorization_.setConstraintStatus(contact_sequence);
  riccati_recursion_.forwardRiccatiRecursionParallel(riccati_factorizer_, 
                                                     ocp.discrete(),
                                                     kkt_matrix, kkt_residual,
                                                     constraint_factorization_);
  const bool exist_state_constraint = ocp.discrete().existStateConstraint();
  if (exist_state_constraint) {
    riccati_recursion_.forwardStateConstraintFactorization(
        riccati_factorizer_, ocp.discrete(), kkt_matrix, kkt_residual, 
        riccati_factorization_);
    riccati_recursion_.backwardStateConstraintFactorization(
        riccati_factorizer_, ocp.discrete(), kkt_matrix, 
        constraint_factorization_);
  }
  direction_calculator_.computeInitialStateDirection<UseContinuationMethod>(
      robots, q, v, s, d, sampling_period);
  if (exist_state_constraint) {
    constraint_factorizer_.computeLagrangeMultiplierDirection(
        ocp.discrete(), riccati_factorization_, constraint_factorization_, d);
    constraint_factorizer_.aggregateLagrangeMultiplierDirection(
       constraint_factorization_,  ocp.discrete(), d, riccati_factorization_);
  }
  riccati_recursion_.forwardRiccatiRecursion(
      riccati_factorizer_, ocp.discrete(), kkt_matrix, kkt_residual, 
      riccati_factorization_, d);
  direction_calculator_.computeNewtonDirectionFromRiccatiFactorization(
      ocp, robots, riccati_factorizer_, riccati_factorization_, s, d);
}

} // namespace idocp 

#endif // IDOCP_RICCATI_SOLVER_HXX_ 