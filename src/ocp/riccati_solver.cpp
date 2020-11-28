#include "idocp/ocp/riccati_solver.hpp"


namespace idocp {

RiccatiSolver::RiccatiSolver(const Robot& robot, const double T, const int N, 
                             const int max_num_impulse, const int num_proc) 
  : riccati_recursion_(robot, T, N, max_num_impulse, num_proc),
    riccati_factorizer_(N, max_num_impulse, robot),
    riccati_factorization_(N, max_num_impulse, robot),
    constraint_factorizer_(robot, N, max_num_impulse, num_proc),
    constraint_factorization_(robot, N, max_num_impulse),
    ocp_direction_calculator_(T, N, max_num_impulse, num_proc) {
}


RiccatiSolver::RiccatiSolver() 
  : riccati_recursion_(),
    riccati_factorizer_(),
    riccati_factorization_(),
    constraint_factorizer_(),
    constraint_factorization_(),
    ocp_direction_calculator_() {
}


RiccatiSolver::~RiccatiSolver() {
}


void RiccatiSolver::setConstraintDimensions(
    const ContactSequence& contact_sequence) {
  constraint_factorization_.setConstraintStatus(contact_sequence);
}


void RiccatiSolver::computeNewtonDirection(
    HybridOCP& split_ocps, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const HybridSolution& s, HybridDirection& d, 
    HybridKKTMatrix& kkt_matrix, HybridKKTResidual& kkt_residual) {
  riccati_recursion_.backwardRiccatiRecursionTerminal(kkt_matrix, kkt_residual, 
                                                      riccati_factorization_);
  riccati_recursion_.backwardRiccatiRecursion(riccati_factorizer_, 
                                              contact_sequence, kkt_matrix, 
                                              kkt_residual, 
                                              riccati_factorization_);
  riccati_recursion_.forwardRiccatiRecursionParallel(riccati_factorizer_, 
                                                     contact_sequence, 
                                                     kkt_matrix, kkt_residual,
                                                     constraint_factorization_);
  if (contact_sequence.existImpulseStage()) {
    riccati_recursion_.forwardStateConstraintFactorization(
        riccati_factorizer_, contact_sequence, kkt_matrix, kkt_residual, 
        riccati_factorization_);
    riccati_recursion_.backwardStateConstraintFactorization(
        riccati_factorizer_, contact_sequence, kkt_matrix, 
        constraint_factorization_);
  }
  ocp_direction_calculator_.computeInitialStateDirection(robots, q, v, s, d);
  if (contact_sequence.existImpulseStage()) {
    constraint_factorizer_.computeLagrangeMultiplierDirection(
        contact_sequence, riccati_factorization_, constraint_factorization_, d);
    constraint_factorizer_.aggregateLagrangeMultiplierDirection(
        constraint_factorization_, contact_sequence, d, riccati_factorization_);
  }
  riccati_recursion_.forwardRiccatiRecursion(
      riccati_factorizer_, contact_sequence, kkt_matrix, kkt_residual, 
      riccati_factorization_, d);
  ocp_direction_calculator_.computeDirection(
      split_ocps, robots, contact_sequence, riccati_factorizer_, 
      riccati_factorization_, s, d);
}


double RiccatiSolver::maxPrimalStepSize(
    const ContactSequence& contact_sequence) const {
  return ocp_direction_calculator_.maxPrimalStepSize(contact_sequence);
}


double RiccatiSolver::maxDualStepSize(
    const ContactSequence& contact_sequence) const {
  return ocp_direction_calculator_.maxDualStepSize(contact_sequence);
}

} // namespace idocp