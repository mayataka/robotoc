#include "robotoc/solver/unconstr_parnmpc_solver.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrParNMPCSolver::UnconstrParNMPCSolver(const UnconstrParNMPC& parnmpc, 
                                             const SolverOptions& solver_options, 
                                             const int nthreads)
  : robots_(nthreads, parnmpc.robot()),
    parnmpc_(parnmpc),
    backward_correction_(parnmpc, nthreads),
    line_search_(parnmpc, nthreads),
    kkt_matrix_(parnmpc.robot(), parnmpc.N()),
    kkt_residual_(parnmpc.robot(), parnmpc.N()),
    s_(parnmpc.robot(), parnmpc.N()),
    d_(parnmpc.robot(), parnmpc.N()),
    N_(parnmpc.N()),
    nthreads_(nthreads),
    T_(parnmpc.T()),
    dt_(parnmpc.T()/parnmpc.N()),
    solver_options_(solver_options),
    solver_statistics_() {
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


UnconstrParNMPCSolver::UnconstrParNMPCSolver() {
}


UnconstrParNMPCSolver::~UnconstrParNMPCSolver() {
}


void UnconstrParNMPCSolver::setSolverOptions(const SolverOptions& solver_options) {
  solver_options_ = solver_options;
}


void UnconstrParNMPCSolver::initConstraints() {
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i < N_-1) {
      parnmpc_[i].initConstraints(robots_[omp_get_thread_num()], i+1, s_[i]);
    }
    else {
      parnmpc_.terminal.initConstraints(robots_[omp_get_thread_num()], i+1, 
                                        s_[i]);
    }
  }
}


void UnconstrParNMPCSolver::initBackwardCorrection(const double t) {
  backward_correction_.initAuxMat(robots_, parnmpc_, t, s_, 
                                  kkt_matrix_, kkt_residual_);
}


void UnconstrParNMPCSolver::updateSolution(const double t, 
                                           const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  backward_correction_.coarseUpdate(robots_, parnmpc_, t, q, v, kkt_matrix_,
                                    kkt_residual_, s_);
  backward_correction_.backwardCorrection(parnmpc_, s_, kkt_matrix_, 
                                          kkt_residual_, d_);
  double primal_step_size     = backward_correction_.primalStepSize();
  const double dual_step_size = backward_correction_.dualStepSize();
  if (solver_options_.enable_line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(parnmpc_, robots_, t, q, v, 
                                                    s_, d_, max_primal_step_size);
  }
  solver_statistics_.primal_step_size.push_back(primal_step_size);
  solver_statistics_.dual_step_size.push_back(dual_step_size);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i < N_-1) {
      parnmpc_[i].updatePrimal(robots_[omp_get_thread_num()], primal_step_size, 
                               d_[i], s_[i]);
      parnmpc_[i].updateDual(dual_step_size);
    }
    else {
      parnmpc_.terminal.updatePrimal(robots_[omp_get_thread_num()],  
                                     primal_step_size, d_[i], s_[i]);
      parnmpc_.terminal.updateDual(dual_step_size);
    }
  }
} 


void UnconstrParNMPCSolver::solve(const double t, const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& v,
                                  const bool init_solver) {
  if (init_solver) {
    initConstraints();
    initBackwardCorrection(t);
    line_search_.clearFilter();
  }
  solver_statistics_.clear(); 
  for (int iter=0; iter<solver_options_.max_iter; ++iter) {
    updateSolution(t, q, v);
    const double kkt_error = KKTError();
    solver_statistics_.kkt_error.push_back(kkt_error); 
    if (kkt_error < solver_options_.kkt_tol) {
      solver_statistics_.convergence = true;
      solver_statistics_.iter = iter+1;
      break;
    }
  }
  if (!solver_statistics_.convergence) {
    solver_statistics_.iter = solver_options_.max_iter;
  }
}


const SolverStatistics& UnconstrParNMPCSolver::getSolverStatistics() const {
  return solver_statistics_;
}


const SplitSolution& UnconstrParNMPCSolver::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= N_);
  return s_[stage];
}


std::vector<Eigen::VectorXd> UnconstrParNMPCSolver::getSolution(
    const std::string& name) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<N_; ++i) {
      sol.push_back(s_[i].q);
    }
  }
  if (name == "v") {
    for (int i=0; i<N_; ++i) {
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


void UnconstrParNMPCSolver::setSolution(const std::string& name, 
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


double UnconstrParNMPCSolver::KKTError(const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  parnmpc_.discretize(t);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i == 0) {
      parnmpc_[0].computeKKTResidual(robots_[omp_get_thread_num()], 
                                     parnmpc_.gridInfo(0), q, v, s_[0], s_[1], 
                                     kkt_matrix_[0], kkt_residual_[0]);
    }
    else if (i < N_-1) {
      parnmpc_[i].computeKKTResidual(robots_[omp_get_thread_num()], 
                                     parnmpc_.gridInfo(i), s_[i-1].q, s_[i-1].v, 
                                     s_[i], s_[i+1], kkt_matrix_[i], kkt_residual_[i]);
    }
    else {
      parnmpc_.terminal.computeKKTResidual(robots_[omp_get_thread_num()],  
                                           parnmpc_.gridInfo(i), s_[i-1].q, 
                                           s_[i-1].v, s_[i], kkt_matrix_[i], 
                                           kkt_residual_[i]);
    }
  }
  return KKTError();
}


double UnconstrParNMPCSolver::KKTError() const {
  double kkt_error = 0;
  for (int i=0; i<N_; ++i) {
    kkt_error += kkt_residual_[i].kkt_error;
  }
  return std::sqrt(kkt_error);
}


double UnconstrParNMPCSolver::cost() const {
  double total_cost = 0;
  for (int i=0; i<N_-1; ++i) {
    total_cost += parnmpc_[i].stageCost();
  }
  total_cost += parnmpc_.terminal.stageCost();
  return total_cost;
}


bool UnconstrParNMPCSolver::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    if (i < N_-1) {
      const bool feasible = parnmpc_[i].isFeasible(robots_[0], s_[i]);
      if (!feasible) {
        std::cout << "INFEASIBLE at time stage " << i << std::endl;
        return false;
      }
    }
    else {
      const bool feasible = parnmpc_.terminal.isFeasible(robots_[0], s_[i]);
      if (!feasible) {
        std::cout << "INFEASIBLE at time stage " << i << std::endl;
        return false;
      }
    }
  }
  return true;
}

} // namespace robotoc