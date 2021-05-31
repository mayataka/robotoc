#include "idocp/solver/unconstr_parnmpc_solver.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

UnconstrParNMPCSolver::UnconstrParNMPCSolver(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints, 
    const double T, const int N, const int nthreads)
  : robots_(nthreads, robot),
    parnmpc_(robot, cost, constraints, N),
    backward_correction_(robot, T, N, nthreads),
    line_search_(robot, T, N, nthreads),
    kkt_matrix_(robot, N),
    kkt_residual_(robot, N),
    s_(robot, N),
    d_(robot, N),
    N_(N),
    nthreads_(nthreads),
    T_(T),
    dt_(T/N),
    kkt_error_(Eigen::VectorXd::Zero(N)) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
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
                                           const Eigen::VectorXd& v, 
                                           const bool line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  backward_correction_.coarseUpdate(robots_, parnmpc_, t, q, v, kkt_matrix_,
                                    kkt_residual_, s_);
  backward_correction_.backwardCorrection(robots_, parnmpc_, s_, 
                                          kkt_matrix_, kkt_residual_, d_);
  double primal_step_size     = backward_correction_.primalStepSize();
  const double dual_step_size = backward_correction_.dualStepSize();
  if (line_search) {
    const double max_primal_step_size = primal_step_size;
    primal_step_size = line_search_.computeStepSize(parnmpc_, robots_, t, q, v, 
                                                    s_, d_, max_primal_step_size);
  }
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


void UnconstrParNMPCSolver::getStateFeedbackGain(const int time_stage, 
                                                 Eigen::MatrixXd& Kq, 
                                                 Eigen::MatrixXd& Kv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  // riccati_solver_.getStateFeedbackGain(time_stage, Kq, Kv);
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


void UnconstrParNMPCSolver::clearLineSearchFilter() {
  line_search_.clearFilter();
}


void UnconstrParNMPCSolver::computeKKTResidual(const double t, 
                                               const Eigen::VectorXd& q, 
                                               const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i == 0) {
      parnmpc_[0].computeKKTResidual(robots_[omp_get_thread_num()], t+dt_,
                                     dt_, q, v, s_[0], s_[1], 
                                     kkt_matrix_[0], kkt_residual_[0]);
    }
    else if (i < N_-1) {
      parnmpc_[i].computeKKTResidual(robots_[omp_get_thread_num()], t+(i+1)*dt_,   
                                     dt_, s_[i-1].q, s_[i-1].v, s_[i], s_[i+1],
                                     kkt_matrix_[i], kkt_residual_[i]);
    }
    else {
      parnmpc_.terminal.computeKKTResidual(robots_[omp_get_thread_num()], t+T_, 
                                           dt_, s_[i-1].q, s_[i-1].v, s_[i],
                                           kkt_matrix_[i], kkt_residual_[i]);
    }
  }
}


double UnconstrParNMPCSolver::KKTError() {
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i < N_-1) {
      kkt_error_.coeffRef(i) 
          = parnmpc_[i].squaredNormKKTResidual(kkt_residual_[0], dt_);
    }
    else {
      kkt_error_.coeffRef(i) 
          = parnmpc_.terminal.squaredNormKKTResidual(kkt_residual_[i], dt_);
    }
  }
  return std::sqrt(kkt_error_.sum());
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

} // namespace idocp