#include "idocp/ocp/riccati_solver.hpp"

#include <stdexcept>


namespace idocp {

RiccatiSolver::RiccatiSolver(const Robot& robot, const int N, 
                             const int max_num_impulse, const int nthreads) 
  : riccati_recursion_(robot, N, max_num_impulse),
    factorization_(robot, N, max_num_impulse),
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
    factorization_(),
    direction_calculator_() {
}


RiccatiSolver::~RiccatiSolver() {
}


void RiccatiSolver::computeNewtonDirection(
    OCP& ocp, std::vector<Robot>& robots, const Eigen::VectorXd& q, 
    const Eigen::VectorXd& v, const Solution& s, Direction& d, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, 
    const StateConstraintJacobian& jac) {
  riccati_recursion_.backwardRiccatiRecursion(ocp.discrete(), kkt_matrix, 
                                              kkt_residual, jac, factorization_);
  direction_calculator_.computeInitialStateDirection(robots, q, v, kkt_matrix, 
                                                     s, d);
  riccati_recursion_.forwardRiccatiRecursion(ocp.discrete(), kkt_matrix, 
                                             kkt_residual, d);
  direction_calculator_.computeNewtonDirection(ocp, robots, factorization_, s, d);
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
}

} // namespace idocp