#include "idocp/ocp/riccati_solver.hpp"

#include <stdexcept>


namespace idocp {

RiccatiSolver::RiccatiSolver(const Robot& robot, const int N, 
                             const int max_num_impulse, const int nthreads) 
  : riccati_recursion_(robot, N, nthreads),
    riccati_factorizer_(robot, N, max_num_impulse),
    riccati_factorization_(robot, N, max_num_impulse),
    constraint_factorizer_(robot, N, max_num_impulse, nthreads),
    constraint_factorization_(robot, N, max_num_impulse),
    direction_calculator_(N, max_num_impulse, nthreads) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


RiccatiSolver::RiccatiSolver() 
  : riccati_recursion_(),
    riccati_factorizer_(),
    riccati_factorization_(),
    constraint_factorizer_(),
    constraint_factorization_(),
    direction_calculator_() {
}


RiccatiSolver::~RiccatiSolver() {
}


void RiccatiSolver::computeNewtonDirection(
    OCP& ocp, std::vector<Robot>& robots, 
    const ContactSequence& contact_sequence, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, Direction& d, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
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
  direction_calculator_.computeInitialStateDirection(robots, q, v, kkt_matrix, 
                                                     s, d);
  riccati_recursion_.forwardRiccatiRecursion(
      riccati_factorizer_, ocp.discrete(), kkt_matrix, kkt_residual, 
      riccati_factorization_, d);
  direction_calculator_.computeNewtonDirectionFromRiccatiFactorization(
      ocp, robots, riccati_factorizer_, riccati_factorization_, s, d);
}


double RiccatiSolver::maxPrimalStepSize() const {
  return direction_calculator_.maxPrimalStepSize();
}


double RiccatiSolver::maxDualStepSize() const {
  return direction_calculator_.maxDualStepSize();
}


void RiccatiSolver::getStateFeedbackGain(const int time_stage, 
                                         Eigen::MatrixXd& Kq, 
                                         Eigen::MatrixXd& Kv) const {
  riccati_factorizer_[time_stage].getStateFeedbackGain(Kq, Kv);
}

} // namespace idocp