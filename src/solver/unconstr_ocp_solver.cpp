#include "robotoc/solver/unconstr_ocp_solver.hpp"

#include <omp.h>
#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrOCPSolver::UnconstrOCPSolver(const UnconstrOCP& ocp, 
                                     const SolverOptions& solver_options, 
                                     const int nthreads)
  : robots_(nthreads, ocp.robot()),
    ocp_(ocp),
    riccati_recursion_(ocp),
    line_search_(ocp, nthreads),
    kkt_matrix_(ocp.N()+1, SplitKKTMatrix(ocp.robot())),
    kkt_residual_(ocp.N()+1, SplitKKTResidual(ocp.robot())),
    s_(ocp.N()+1, SplitSolution(ocp.robot())),
    d_(ocp.N()+1, SplitDirection(ocp.robot())),
    riccati_factorization_(ocp.N()+1, SplitRiccatiFactorization(ocp.robot())),
    N_(ocp.N()),
    nthreads_(nthreads),
    T_(ocp.T()),
    dt_(ocp.T()/ocp.N()),
    primal_step_size_(Eigen::VectorXd::Zero(ocp.N())), 
    dual_step_size_(Eigen::VectorXd::Zero(ocp.N())),
    solver_options_(solver_options),
    solver_statistics_() {
  if (nthreads <= 0) {
    throw std::out_of_range("[UnconstrOCPSolver] invalid argument: nthreads must be positive!");
  }
  initConstraints();
}


UnconstrOCPSolver::UnconstrOCPSolver() {
}


UnconstrOCPSolver::~UnconstrOCPSolver() {
}


void UnconstrOCPSolver::setSolverOptions(const SolverOptions& solver_options) {
  solver_options_ = solver_options;
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
                                       const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  ocp_.discretize(t);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i == 0) {
      ocp_[0].computeKKTSystem(robots_[omp_get_thread_num()], ocp_.gridInfo(0),  
                               s_[0], s_[1], kkt_matrix_[0], kkt_residual_[0]);
    }
    else if (i < N_) {
      ocp_[i].computeKKTSystem(robots_[omp_get_thread_num()], ocp_.gridInfo(i),
                               s_[i], s_[i+1], kkt_matrix_[i], kkt_residual_[i]);
    }
    else {
      ocp_.terminal.computeKKTSystem(robots_[omp_get_thread_num()], 
                                     ocp_.gridInfo(N_), s_[N_-1].q, s_[N_], 
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
  if (solver_options_.enable_line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(ocp_, robots_, t, q, v, s_,
                                                    d_, max_primal_step_size);
  }
  solver_statistics_.primal_step_size.push_back(primal_step_size);
  solver_statistics_.dual_step_size.push_back(dual_step_size);
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


void UnconstrOCPSolver::solve(const double t, const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v,
                              const bool init_solver) {
  if (q.size() != robots_[0].dimq()) {
    throw std::out_of_range("[UnconstrOCPSolver] invalid argument: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
  }
  if (v.size() != robots_[0].dimv()) {
    throw std::out_of_range("[UnconstrOCPSolver] invalid argument: v.size() must be " + std::to_string(robots_[0].dimq()) + "!");
  }
  if (solver_options_.enable_benchmark) {
    timer_.tick();
  }
  if (init_solver) {
    initConstraints();
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
  if (solver_options_.enable_benchmark) {
    timer_.tock();
    solver_statistics_.cpu_time = timer_.ms();
  }
}


const SolverStatistics& UnconstrOCPSolver::getSolverStatistics() const {
  return solver_statistics_;
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


const std::vector<LQRPolicy>& UnconstrOCPSolver::getLQRPolicy() const {
  return riccati_recursion_.getLQRPolicy();
}


void UnconstrOCPSolver::setSolution(const std::string& name, 
                                    const Eigen::VectorXd& value) {
  if (name == "q") {
    for (auto& e : s_) { e.q = value; }
  }
  else if (name == "v") {
    for (auto& e : s_) { e.v = value; }
  }
  else if (name == "a") {
    for (auto& e : s_) { e.a  = value; }
  }
  else if (name == "u") {
    for (auto& e : s_) { e.u = value; }
  }
  else {
    throw std::invalid_argument("[UnconstrOCPSolver] invalid arugment: name must be q, v, a, or u!");
  }
  initConstraints();
}


double UnconstrOCPSolver::KKTError(const double t, const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v) {
  if (q.size() != robots_[0].dimq()) {
    throw std::out_of_range("[UnconstrOCPSolver] invalid argument: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
  }
  if (v.size() != robots_[0].dimv()) {
    throw std::out_of_range("[UnconstrOCPSolver] invalid argument: v.size() must be " + std::to_string(robots_[0].dimv()) + "!");
  }
  ocp_.discretize(t);
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i == 0) {
      ocp_[0].computeKKTResidual(robots_[omp_get_thread_num()], ocp_.gridInfo(0), 
                                 s_[0], s_[1], kkt_matrix_[0], kkt_residual_[0]);
    }
    else if (i < N_) {
      ocp_[i].computeKKTResidual(robots_[omp_get_thread_num()], ocp_.gridInfo(i),   
                                 s_[i], s_[i+1], kkt_matrix_[i], kkt_residual_[i]);
    }
    else {
      ocp_.terminal.computeKKTResidual(robots_[omp_get_thread_num()],  
                                       ocp_.gridInfo(N_), s_[N_-1].q, s_[N_], 
                                       kkt_matrix_[N_], kkt_residual_[N_]);
    }
  }
  return KKTError();
}


double UnconstrOCPSolver::KKTError() const {
  double kkt_error = 0;
  for (int i=0; i<=N_; ++i) {
    kkt_error += kkt_residual_[i].KKTError();
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


void UnconstrOCPSolver::setRobotProperties(const RobotProperties& properties) {
  for (auto& e : robots_) {
    e.setRobotProperties(properties);
  }
}

} // namespace robotoc