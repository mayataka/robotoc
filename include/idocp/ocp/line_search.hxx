#include "idocp/ocp/line_search.hpp"

template <typename OCPType>
inline double computeStepSize(const std::vector<SplitOCPType>& split_ocp,  
                              const double t, const Eigen::VectorXd& q, 
                              const Eigen::VectorXd& v, 
                              const std::vector<SplitSolution>& s,
                              const std::vector<SplitDirection>& d,
                              const int max_iteration) {
  // If filter is empty, augment the current solution to the filter.
  if (filter_.isEmpty()) {
    #pragma omp parallel for num_threads(num_proc_)
    for (int i=0; i<=N_; ++i) {
      if (i < N_) {
        const int robot_id = omp_get_thread_num();
        costs_.coeffRef(i) = split_ocps_[i].stageCost(robots_[robot_id], t+i*dtau_, dtau_, s_[i]);
        violations_.coeffRef(i) = split_ocps_[i].constraintViolation(
            robots_[robot_id], contact_sequence_.contactStatus(i), t+i*dtau_, 
            dtau_, s_[i], s_[i+1].q, s_[i+1].v);
      }
      else {
        const int robot_id = omp_get_thread_num();
        costs_.coeffRef(N_) = terminal_ocp_.terminalCost(robots_[robot_id], t+T_, s_[N_]);
      }
    }
    filter_.augment(costs_.sum(), violations_.sum());
  }
  while (primal_step_size > min_step_size_) {
    #pragma omp parallel for num_threads(num_proc_)
    for (int i=0; i<=N_; ++i) {
      if (i < N_) {
        const int robot_id = omp_get_thread_num();
        s_tmp_[i].setTemporarySolution(robots_[robot_id], contact_sequence_.contactStatus(i), 
                                        primal_step_size, s_[i], d_[i], s_[i+1], d_[i+1]);
        costs_.coeffRef(i) = split_ocps_[i].stageCost(robots_[robot_id], t+i*dtau_, dtau_, 
                                                      s_tmp_[i].splitSolution(), 
                                                      primal_step_size);
        violations_.coeffRef(i) = split_ocps_[i].constraintViolation(
            robots_[robot_id], contact_sequence_.contactStatus(i), t+i*dtau_, 
            dtau_, s_tmp_[i].splitSolution(), s_tmp_[i].q_next(), s_tmp_[i].v_next());
      }
      else {
        const int robot_id = omp_get_thread_num();
        s_tmp_[N_].setTemporarySolution(robots_[robot_id], primal_step_size, s_[N_], d_[N_]);
        costs_.coeffRef(N_) = terminal_ocp_.terminalCost(robots_[robot_id], t+T_, s_tmp_[N_].splitSolution());
      }
    }
    const double cost_sum = costs_.sum();
    const double violation_sum = violations_.sum();
    if (filter_.isAccepted(cost_sum, violation_sum)) {
      filter_.augment(cost_sum, violation_sum);
      break;
    }
    primal_step_size *= step_size_reduction_rate_;
  }
}