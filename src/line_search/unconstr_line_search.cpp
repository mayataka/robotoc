#include "idocp/line_search/unconstr_line_search.hpp"

#include <stdexcept>
#include <cassert>


namespace idocp {

UnconstrLineSearch::UnconstrLineSearch(const Robot& robot, const double T, 
                                       const int N, const int nthreads, 
                                       const double step_size_reduction_rate,
                                       const double min_step_size) 
  : filter_(),
    N_(N), 
    nthreads_(nthreads),
    T_(T),
    dt_(T/N),
    step_size_reduction_rate_(step_size_reduction_rate), 
    min_step_size_(min_step_size),
    costs_(Eigen::VectorXd::Zero(N+1)), 
    violations_(Eigen::VectorXd::Zero(N)), 
    s_try_(robot, N), 
    kkt_residual_(robot, N) {
}


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
    s_try_(), 
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
      costs_.coeffRef(i) = ocp[i].stageCost(robots[omp_get_thread_num()], 
                                            t+i*dt_, dt_, s[i], 
                                            primal_step_size);
      violations_.coeffRef(i) = ocp[i].constraintViolation(
          robots[omp_get_thread_num()], t+i*dt_, dt_, s[i], s[i+1].q, s[i+1].v,
          kkt_residual_[i]);
    }
    else {
      costs_.coeffRef(i) = ocp.terminal.terminalCost(robots[omp_get_thread_num()], 
                                                     t+T_, s[i]);
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
      costs_.coeffRef(0) = parnmpc[0].stageCost(robots[omp_get_thread_num()], 
                                                t+dt_, dt_, s[0], 
                                                primal_step_size);
      violations_.coeffRef(0) = parnmpc[0].constraintViolation(
          robots[omp_get_thread_num()], t+dt_, dt_, q, v, s[0], 
          kkt_residual_[0]);
    }
    else if (i < N_-1) {
      costs_.coeffRef(i) = parnmpc[i].stageCost(robots[omp_get_thread_num()], 
                                                t+(i+1)*dt_, dt_, s[i], 
                                                primal_step_size);
      violations_.coeffRef(i) = parnmpc[i].constraintViolation(
          robots[omp_get_thread_num()], t+(i+1)*dt_, dt_, 
          s[i-1].q, s[i-1].v, s[i], kkt_residual_[i]);
    }
    else {
      costs_.coeffRef(i) = parnmpc.terminal.stageCost(
          robots[omp_get_thread_num()], t+(i+1)*dt_, dt_, s[i], 
          primal_step_size);
      violations_.coeffRef(i) = parnmpc.terminal.constraintViolation(
          robots[omp_get_thread_num()], t+(i+1)*dt_, dt_, s[i-1].q, 
          s[i-1].v, s[i], kkt_residual_[i]);
    }
  }
}


void UnconstrLineSearch::computeTrySolution(const Solution& s, 
                                            const Direction& d, 
                                            const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  for (int i=0; i<=N_; ++i) {
    computeTrySolution(s[i], d[i], step_size, s_try_[i]);
  }
}

} // namespace idocp