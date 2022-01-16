#include "robotoc/line_search/unconstr_line_search.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

UnconstrLineSearch::UnconstrLineSearch() 
  : filter_(),
    N_(0), 
    nthreads_(0),
    T_(0),
    dt_(0),
    step_size_reduction_rate_(0), 
    min_step_size_(0),
    costs_(), 
    violations_(), 
    s_trial_(), 
    kkt_residual_() {
}


UnconstrLineSearch::~UnconstrLineSearch() {
}


void UnconstrLineSearch::clearFilter() {
  filter_.clear();
}


bool UnconstrLineSearch::isFilterEmpty() const {
  return filter_.isEmpty();
}


void UnconstrLineSearch::computeCostAndViolation(UnconstrOCP& ocp, 
                                                 aligned_vector<Robot>& robots, 
                                                 const double t, 
                                                 const Eigen::VectorXd& q, 
                                                 const Eigen::VectorXd& v, 
                                                 const Solution& s, 
                                                 const double primal_step_size) {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  clearCosts();
  clearViolations();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      ocp[i].evalOCP(robots[omp_get_thread_num()], i, t+i*dt_, dt_, s[i], 
                     s[i+1].q, s[i+1].v, kkt_residual_[i]);
      costs_.coeffRef(i) = ocp[i].stageCost();
      violations_.coeffRef(i) = ocp[i].constraintViolation(kkt_residual_[i], dt_);
    }
    else {
      ocp.terminal.evalOCP(robots[omp_get_thread_num()], t+i*dt_, 
                           s[i], kkt_residual_[i]);
      costs_.coeffRef(i) = ocp.terminal.terminalCost();
    }
  }
}


void UnconstrLineSearch::computeCostAndViolation(UnconstrParNMPC& parnmpc, 
                                                 aligned_vector<Robot>& robots, 
                                                 const double t, 
                                                 const Eigen::VectorXd& q, 
                                                 const Eigen::VectorXd& v, 
                                                 const Solution& s, 
                                                 const double primal_step_size) {
  assert(robots.size() == nthreads_);
  assert(q.size() == robots[0].dimq());
  assert(v.size() == robots[0].dimv());
  clearCosts();
  clearViolations();
  #pragma omp parallel for num_threads(nthreads_)
  for (int i=0; i<N_; ++i) {
    if (i == 0) {
      parnmpc[0].evalOCP(robots[omp_get_thread_num()], 0, t+dt_, dt_, q, v, 
                         s[0], kkt_residual_[0]);
      costs_.coeffRef(0) = parnmpc[0].stageCost();
      violations_.coeffRef(0) = parnmpc[0].constraintViolation(kkt_residual_[0], dt_);
    }
    else if (i < N_-1) {
      parnmpc[i].evalOCP(robots[omp_get_thread_num()], i, t+(i+1)*dt_, dt_, 
                         s[i-1].q, s[i-1].v, s[i], kkt_residual_[i]);
      costs_.coeffRef(i) = parnmpc[i].stageCost();
      violations_.coeffRef(i) = parnmpc[i].constraintViolation(kkt_residual_[i], dt_);
    }
    else {
      parnmpc.terminal.evalOCP(robots[omp_get_thread_num()], i, t+(i+1)*dt_, dt_, 
                               s[i-1].q, s[i-1].v, s[i], kkt_residual_[i]);
      costs_.coeffRef(i) = parnmpc.terminal.stageCost();
      violations_.coeffRef(i) = parnmpc.terminal.constraintViolation(kkt_residual_[i], dt_);
    }
  }
}


void UnconstrLineSearch::computeSolutionTrial(const Solution& s, 
                                              const Direction& d, 
                                              const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  for (int i=0; i<=N_; ++i) {
    computeSolutionTrial(s[i], d[i], step_size, s_trial_[i]);
  }
}

} // namespace robotoc