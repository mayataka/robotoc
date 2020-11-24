#include "idocp/ocp/ocp_riccati_solver.hpp"


namespace idocp {

OCPRiccatiSolver::OCPRiccatiSolver(const Robot& robot, const double T, 
                                   const int N, 
                                   const int max_num_impulse, 
                                   const int num_proc) 
  : riccati_recursion_(robot, T, N, max_num_impulse, num_proc),
    riccati_factorization_(N, max_num_impulse, robot),
    constraint_factorizer_(robot, max_num_impulse, num_proc),
    constraint_factorization_(
        max_num_impulse, StateConstraintRiccatiFactorization(robot, N, 
                                                             max_num_impulse)),
    ocp_direction_calculator_(T, N, max_num_impulse, num_proc) {
}


OCPRiccatiSolver::OCPRiccatiSolver() 
  : riccati_recursion_(),
    riccati_factorization_(),
    constraint_factorizer_(),
    constraint_factorization_(),
    ocp_direction_calculator_() {
}


OCPRiccatiSolver::~OCPRiccatiSolver() {
}


void OCPRiccatiSolver::setConstraintDimensions(
    const ContactSequence& contact_sequence) {
  const int num_impulse_stages = contact_sequence.totalNumImpulseStages();
  for (int i=0; i<num_impulse_stages; ++i) {
    constraint_factorization_[i].setImpulseStatus(
        contact_sequence.impulseStatus(i));
  }
}


void OCPRiccatiSolver::computeDirection(HybridOCP& split_ocps, 
                                        std::vector<Robot>& robots, 
                                        const ContactSequence& contact_sequence, 
                                        const Eigen::VectorXd& q, 
                                        const Eigen::VectorXd& v, 
                                        const HybridSolution& s, 
                                        HybridDirection& d, 
                                        HybridKKTMatrix& kkt_matrix, 
                                        HybridKKTResidual& kkt_residual) {
  riccati_recursion_.backwardRiccatiRecursionTerminal(kkt_matrix, 
                                                      kkt_residual, 
                                                      riccati_factorization_);
  riccati_recursion_.backwardRiccatiRecursion(contact_sequence, kkt_matrix, 
                                              kkt_residual, 
                                              riccati_factorization_);
  riccati_recursion_.forwardRiccatiRecursionParallel(contact_sequence, 
                                                     kkt_matrix, kkt_residual,
                                                     constraint_factorization_);
  riccati_recursion_.forwardRiccatiRecursionSerial(contact_sequence, 
                                                   kkt_matrix, kkt_residual, 
                                                   riccati_factorization_);
  riccati_recursion_.backwardStateConstraintFactorization(
      contact_sequence, kkt_matrix, constraint_factorization_);
  ocp_direction_calculator_.computeInitialStateDirection(robots, q, v, s, d);
  // if (contact_sequence.existImpulseStage()) {
  //   constraint_factorizer_.computeLagrangeMultiplierDirection(
  //       contact_sequence, riccati_factorization_.impulse, 
  //       constraint_factorization_, d[0].dx(), d.impulse);
  //   ocp_direction_calculator_.aggregateLagrangeMultiplierDirection(
  //       contact_sequence, riccati_factorization_, constraint_factorization_, d);
  // }
  ocp_direction_calculator_.computeDirection(
      split_ocps, robots, contact_sequence, 
      riccati_recursion_.getFactorizersHandle(), riccati_factorization_, 
      constraint_factorization_, s, d);
}


double OCPRiccatiSolver::maxPrimalStepSize(
    const ContactSequence& contact_sequence) const {
  return ocp_direction_calculator_.maxPrimalStepSize(contact_sequence);
}


double OCPRiccatiSolver::maxDualStepSize(
    const ContactSequence& contact_sequence) const {
  return ocp_direction_calculator_.maxDualStepSize(contact_sequence);
}

} // namespace idocp