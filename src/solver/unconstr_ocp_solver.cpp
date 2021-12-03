#include "robotoc/solver/unconstr_ocp_solver.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrOCPSolver::UnconstrOCPSolver(const UnconstrOCP& ocp, const int nthreads)
  : robots_(nthreads, ocp.robot()),
    ocp_(ocp),
    riccati_recursion_(ocp),
    line_search_(ocp.robot(), ocp.T(), ocp.N(), nthreads),
    kkt_matrix_(ocp.robot(), ocp.N()),
    kkt_residual_(ocp.robot(), ocp.N()),
    s_(ocp.robot(), ocp.N()),
    d_(ocp.robot(), ocp.N()),
    riccati_factorization_(ocp.N()+1, SplitRiccatiFactorization(ocp.robot())),
    N_(ocp.N()),
    nthreads_(nthreads),
    T_(ocp.T()),
    dt_(ocp.T()/ocp.N()),
    primal_step_size_(Eigen::VectorXd::Zero(ocp.N())), 
    dual_step_size_(Eigen::VectorXd::Zero(ocp.N())) {
  try {
    if (nthreads <= 0) {
      throw std::out_of_range("invalid value: nthreads must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  initConstraints();
}


UnconstrOCPSolver::UnconstrOCPSolver() {
}


UnconstrOCPSolver::~UnconstrOCPSolver() {
}


void UnconstrOCPSolver::initConstraints() {
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      ocp_[i].initConstraints(robots_[omp_get_thread_num()], i, s_[i]);
    }
    else {
      ocp_.terminal.initConstraints(robots_[omp_get_thread_num()], N_, s_[N_]);
    }
  }
}


void UnconstrOCPSolver::updateSolution(const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v, 
                                       const bool line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i == 0) {
      ocp_[0].computeKKTSystem(robots_[omp_get_thread_num()], t, dt_, s_[0], 
                               s_[1], kkt_matrix_[0], kkt_residual_[0]);
    }
    else if (i < N_) {
      ocp_[i].computeKKTSystem(robots_[omp_get_thread_num()], t+i*dt_, dt_, s_[i], 
                               s_[i+1], kkt_matrix_[i], kkt_residual_[i]);
    }
    else {
      ocp_.terminal.computeKKTSystem(robots_[omp_get_thread_num()], t+T_, 
                                     s_[N_-1].q, s_[N_], 
                                     kkt_matrix_[N_], kkt_residual_[N_]);
    }
  }
  riccati_recursion_.backwardRiccatiRecursion(kkt_matrix_, kkt_residual_,
                                              riccati_factorization_);
  d_[0].dq() = q - s_[0].q;
  d_[0].dv() = v - s_[0].v;
  riccati_recursion_.forwardRiccatiRecursion(kkt_residual_, d_);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    UnconstrRiccatiFactorizer::computeCostateDirection(riccati_factorization_[i], 
                                                       d_[i]);
    if (i < N_) {
      ocp_[i].expandPrimalAndDual(dt_, s_[i], kkt_matrix_[i], 
                                  kkt_residual_[i], d_[i]);
      primal_step_size_.coeffRef(i) = ocp_[i].maxPrimalStepSize();
      dual_step_size_.coeffRef(i)   = ocp_[i].maxDualStepSize();
    }
  }
  double primal_step_size = primal_step_size_.minCoeff();
  const double dual_step_size   = dual_step_size_.minCoeff();
  if (line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(ocp_, robots_, t, q, v, s_,
                                                    d_, max_primal_step_size);
  }
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      ocp_[i].updatePrimal(robots_[omp_get_thread_num()], primal_step_size, 
                           d_[i], s_[i]);
      ocp_[i].updateDual(dual_step_size);
    }
    else {
      ocp_.terminal.updatePrimal(robots_[omp_get_thread_num()],  
                                 primal_step_size, d_[N_], s_[N_]);
      ocp_.terminal.updateDual(dual_step_size);
    }
  }
} 


const SplitSolution& UnconstrOCPSolver::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= N_);
  return s_[stage];
}


std::vector<Eigen::VectorXd> UnconstrOCPSolver::getSolution(
    const std::string& name) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<=N_; ++i) {
      sol.push_back(s_[i].q);
    }
  }
  if (name == "v") {
    for (int i=0; i<=N_; ++i) {
      sol.push_back(s_[i].v);
    }
  }
  if (name == "a") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].a);
    }
  }
  if (name == "u") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].u);
    }
  }
  return sol;
}


void UnconstrOCPSolver::getStateFeedbackGain(const int time_stage, 
                                             Eigen::MatrixXd& Kq, 
                                             Eigen::MatrixXd& Kv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  riccati_recursion_.getStateFeedbackGain(time_stage, Kq, Kv);
}


void UnconstrOCPSolver::setSolution(const std::string& name, 
                                    const Eigen::VectorXd& value) {
  try {
    if (name == "q") {
      for (auto& e : s_.data) { e.q = value; }
    }
    else if (name == "v") {
      for (auto& e : s_.data) { e.v = value; }
    }
    else if (name == "a") {
      for (auto& e : s_.data) { e.a  = value; }
    }
    else if (name == "u") {
      for (auto& e : s_.data) { e.u = value; }
    }
    else {
      throw std::invalid_argument("invalid arugment: name must be q, v, a, or u!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  initConstraints();
}


void UnconstrOCPSolver::clearLineSearchFilter() {
  line_search_.clearFilter();
}


void UnconstrOCPSolver::computeKKTResidual(const double t, 
                                           const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i == 0) {
      ocp_[0].computeKKTResidual(robots_[omp_get_thread_num()], t, dt_, 
                                 s_[0], s_[1], kkt_matrix_[0], kkt_residual_[0]);
    }
    else if (i < N_) {
      ocp_[i].computeKKTResidual(robots_[omp_get_thread_num()], t+i*dt_, dt_,  
                                 s_[i], s_[i+1], kkt_matrix_[i], kkt_residual_[i]);
    }
    else {
      ocp_.terminal.computeKKTResidual(robots_[omp_get_thread_num()], t+T_, 
                                       s_[N_-1].q, s_[N_], 
                                       kkt_matrix_[N_], kkt_residual_[N_]);
    }
  }
}


double UnconstrOCPSolver::KKTError() {
  double kkt_error = 0;
  for (int i=0; i<=N_; ++i) {
    kkt_error += kkt_residual_[i].kkt_error;
  }
  return std::sqrt(kkt_error);
}


double UnconstrOCPSolver::cost() const {
  double total_cost = 0;
  for (int i=0; i<N_; ++i) {
    total_cost += ocp_[i].stageCost();
  }
  total_cost += ocp_.terminal.terminalCost();
  return total_cost;
}


bool UnconstrOCPSolver::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    const bool feasible = ocp_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time stage " << i << std::endl;
      return false;
    }
  }
  return true;
}

} // namespace robotoc