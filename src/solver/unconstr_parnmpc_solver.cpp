#include "robotoc/solver/unconstr_parnmpc_solver.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrParNMPCSolver::UnconstrParNMPCSolver(const OCP& ocp, 
                                             const SolverOptions& solver_options, 
                                             const int nthreads)
  : robots_(nthreads, ocp.robot),
    time_discretization_(ocp.N+1, GridInfo()),
    backward_correction_(ocp, nthreads),
    line_search_(ocp, nthreads),
    ocp_(ocp),
    kkt_matrix_(ocp.N+1, SplitKKTMatrix(ocp.robot)),
    kkt_residual_(ocp.N+1, SplitKKTResidual(ocp.robot)),
    s_(ocp.N+1, SplitSolution(ocp.robot)),
    d_(ocp.N+1, SplitDirection(ocp.robot)),
    solver_options_(solver_options),
    solver_statistics_() {
  if (!ocp.cost) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: ocp.cost should not be nullptr!");
  }
  if (!ocp.constraints) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: ocp.constraints should not be nullptr!");
  }
  if (ocp.T <= 0) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: ocp.T must be positive!");
  }
  if (ocp.N <= 0) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: ocp.N must be positive!");
  }
  if (nthreads <= 0) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: nthreads must be positive!");
  }
  const double dt = ocp.T / ocp.N;
  for (int i=0; i<=ocp.N; ++i) {
    time_discretization_[i].t = dt * i;
    time_discretization_[i].dt = dt;
    time_discretization_[i].phase = -1;
    time_discretization_[i].stage = i;
    time_discretization_[i].impact_index = -1;
    time_discretization_[i].lift_index = -1;
    time_discretization_[i].stage_in_phase = i;
    time_discretization_[i].num_grids_in_phase = ocp.N;
  }
  time_discretization_[ocp.N].type = GridType::Terminal;
  initConstraints();
}


UnconstrParNMPCSolver::UnconstrParNMPCSolver() 
  : robots_(),
    time_discretization_(),
    backward_correction_(),
    line_search_(),
    ocp_(),
    kkt_matrix_(),
    kkt_residual_(),
    s_(),
    d_(),
    solver_options_(),
    solver_statistics_() {
}


void UnconstrParNMPCSolver::setSolverOptions(const SolverOptions& solver_options) {
  solver_options_ = solver_options;
}


void UnconstrParNMPCSolver::discretize(const double t) {
  const double dt = ocp_.T / ocp_.N;
  for (int i=0; i<=ocp_.N; ++i) {
    time_discretization_[i].t = dt * i;
  }
}


void UnconstrParNMPCSolver::initConstraints() {
  backward_correction_.initConstraints(robots_, time_discretization_, s_);
}


void UnconstrParNMPCSolver::initBackwardCorrection() {
  backward_correction_.initAuxMat(robots_, time_discretization_, s_, 
                                  kkt_matrix_, kkt_residual_);
}


void UnconstrParNMPCSolver::updateSolution(const double t, 
                                           const Eigen::VectorXd& q, 
                                           const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  discretize(t);
  backward_correction_.coarseUpdate(robots_, time_discretization_, q, v, 
                                    s_, kkt_matrix_, kkt_residual_);
  backward_correction_.backwardCorrection(time_discretization_, s_, 
                                          kkt_matrix_, kkt_residual_, d_);
  double primal_step_size     = backward_correction_.primalStepSize();
  const double dual_step_size = backward_correction_.dualStepSize();
  if (solver_options_.enable_line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(backward_correction_, 
                                                    robots_, time_discretization_, 
                                                    q, v, s_, d_, max_primal_step_size);
  }
  solver_statistics_.primal_step_size.push_back(primal_step_size);
  solver_statistics_.dual_step_size.push_back(dual_step_size);
  backward_correction_.integrateSolution(robots_, time_discretization_, 
                                         primal_step_size, dual_step_size, d_, s_);
} 


void UnconstrParNMPCSolver::solve(const double t, const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& v,
                                  const bool init_solver) {
  if (q.size() != robots_[0].dimq()) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
  }
  if (v.size() != robots_[0].dimv()) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: v.size() must be " + std::to_string(robots_[0].dimv()) + "!");
  }
  if (solver_options_.enable_benchmark) {
    timer_.tick();
  }
  if (init_solver) {
    discretize(t);
    initConstraints();
    initBackwardCorrection();
    line_search_.clearFilter();
  }
  solver_statistics_.clear(); 
  for (int iter=0; iter<solver_options_.max_iter; ++iter) {
    updateSolution(t, q, v);
    solver_statistics_.performance_index.push_back(backward_correction_.getEval()); 
    const double kkt_error = KKTError();
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


const SolverStatistics& UnconstrParNMPCSolver::getSolverStatistics() const {
  return solver_statistics_;
}


const SplitSolution& UnconstrParNMPCSolver::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= ocp_.N);
  return s_[stage];
}


std::vector<Eigen::VectorXd> UnconstrParNMPCSolver::getSolution(
    const std::string& name) const {
  std::vector<Eigen::VectorXd> sol;
  if (name == "q") {
    for (int i=0; i<ocp_.N; ++i) {
      sol.push_back(s_[i].q);
    }
  }
  if (name == "v") {
    for (int i=0; i<ocp_.N; ++i) {
      sol.push_back(s_[i].v);
    }
  }
  if (name == "a") {
    for (int i=0; i<ocp_.N; ++i) {
      sol.push_back(s_[i].a);
    }
  }
  if (name == "u") {
    for (int i=0; i<ocp_.N; ++i) {
      sol.push_back(s_[i].u);
    }
  }
  return sol;
}


void UnconstrParNMPCSolver::setSolution(const std::string& name, 
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
    throw std::invalid_argument("[UnconstrParNMPCSolver] invalid arugment: name must be q, v, a, or u!");
  }
  initConstraints();
}


double UnconstrParNMPCSolver::KKTError(const double t, const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v) {
  if (q.size() != robots_[0].dimq()) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: q.size() must be " + std::to_string(robots_[0].dimq()) + "!");
  }
  if (v.size() != robots_[0].dimv()) {
    throw std::out_of_range("[UnconstrParNMPCSolver] invalid argument: v.size() must be " + std::to_string(robots_[0].dimv()) + "!");
  }
  discretize(t);
  backward_correction_.evalKKT(robots_, time_discretization_, q, v, 
                               s_, kkt_matrix_, kkt_residual_);
  return KKTError();
}


double UnconstrParNMPCSolver::KKTError() const {
  return std::sqrt(backward_correction_.getEval().kkt_error);
}


const std::vector<GridInfo>& UnconstrParNMPCSolver::getTimeDiscretization() const {
  return time_discretization_;
}


void UnconstrParNMPCSolver::setRobotProperties(const RobotProperties& properties) {
  for (auto& e : robots_) {
    e.setRobotProperties(properties);
  }
}

} // namespace robotoc