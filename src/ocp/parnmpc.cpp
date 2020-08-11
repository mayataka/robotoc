#include "idocp/ocp/parnmpc.hpp"

#include <utility>
#include <cmath>
#include <omp.h>
#include <assert.h>


namespace idocp {

ParNMPC::ParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                 const std::shared_ptr<Constraints>& constraints, 
                 const double T, const int N, const int num_proc)
  : split_ocps_(N, SplitParNMPC(robot, cost, constraints)),
    robots_(num_proc, robot),
    filter_(),
    T_(T),
    dtau_(T/N),
    step_size_reduction_rate_(0.75),
    min_step_size_(0.05),
    N_(N),
    num_proc_(num_proc),
    s_(N+1, SplitSolution(robot)),
    s_new_(N, SplitSolution(robot)),
    s_old_(N, SplitSolution(robot)),
    aux_mat_(N, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    aux_mat_old_(N+1, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    phiq_(Eigen::VectorXd::Zero(robot.dimv())),
    phiv_(Eigen::VectorXd::Zero(robot.dimv())),
    primal_step_sizes_(Eigen::VectorXd::Zero(N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N)),
    costs_(Eigen::VectorXd::Zero(N)), 
    constraints_violations_(Eigen::VectorXd::Zero(N)),
    contact_sequence_(N, std::vector<bool>(robot.max_point_contacts(), false)) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N; ++i) {
    robot.normalizeConfiguration(s_[i].q);
    robot.normalizeConfiguration(s_new_[i].q);
    robot.normalizeConfiguration(s_old_[i].q);
  }
  bool feasible = isCurrentSolutionFeasible();
  initConstraints();
}


ParNMPC::ParNMPC()
  : split_ocps_(),
    robots_(),
    filter_(),
    T_(0),
    dtau_(0),
    step_size_reduction_rate_(0),
    min_step_size_(0),
    N_(0),
    num_proc_(0),
    s_(),
    s_new_(),
    s_old_(),
    aux_mat_(),
    aux_mat_old_(),
    phiq_(),
    phiv_(),
    primal_step_sizes_(),
    dual_step_sizes_(),
    costs_(), 
    constraints_violations_(),
    contact_sequence_() {
}


ParNMPC::~ParNMPC() {
}


void ParNMPC::updateSolution(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, 
                             const bool use_line_search) {
  int time_step;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    robots_[robot_id].setContactStatus(contact_sequence_[i]);
    bool is_terminal_ocp;
    if (time_step == 0) {
      split_ocps_[i].coarseUpdate(robots_[robot_id], t+(i+1)*dtau_, dtau_, 
                                  q, v, s_[i], s_[i+1].lmd, s_[i+1].gmm, 
                                  s_[i+1].q, aux_mat_old_[i+1], aux_mat_[i], 
                                  s_new_[i], false);
    }
    else if (i < N_-1) {
      split_ocps_[i].coarseUpdate(robots_[robot_id], t+(i+1)*dtau_, dtau_, 
                                  s_[i-1].q, s_[i-1].v, s_[i], 
                                  s_[i+1].lmd, s_[i+1].gmm, s_[i+1].q, 
                                  aux_mat_old_[i+1], aux_mat_[time_step], 
                                  s_new_[time_step], false);
    }
    else {
      split_ocps_[i].coarseUpdate(robots_[robot_id], t+(i+1)*dtau_, dtau_, 
                                  s_[i-1].q, s_[i-1].v, s_[i], 
                                  s_[i+1].lmd, s_[i+1].gmm, s_[i+1].q, 
                                  aux_mat_old_[i+1], aux_mat_[time_step], 
                                  s_new_[time_step], true);
    }
  }
  for (int i=N_-2; i>=0; --i) {
    split_ocps_[i].backwardCollectionSerial(s_old_[i+1], s_new_[i+1], s_new_[i]);
  }
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    robots_[robot_id].setContactStatus(contact_sequence_[i]);
    split_ocps_[i].backwardCollectionParallel(robots_[robot_id], s_new_[i]);
  }
  for (int i=1; i<N_; ++i) {
    split_ocps_[i].forwardCollectionSerial(robots_[0], s_old_[i-1], s_new_[i-1], 
                                           s_new_[i]);
  }
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    robots_[robot_id].setContactStatus(contact_sequence_[i]);
    split_ocps_[i].forwardCollectionParallel(robots_[robot_id], s_new_[i]);
    split_ocps_[i].computePrimalAndDualDirection(robots_[robot_id], dtau_, 
                                                  s_[i], s_new_[i]);
    primal_step_sizes_.coeffRef(i) = split_ocps_[i].maxPrimalStepSize();
    dual_step_sizes_.coeffRef(i) = split_ocps_[i].maxDualStepSize();
  }
  double primal_step_size = primal_step_sizes_.minCoeff();
  const double dual_step_size = dual_step_sizes_.minCoeff();
  // if (use_line_search) {
  //   // If filter is empty, augment the current solution to the filter.
  //   if (filter_.isEmpty()) {
  //     #pragma omp parallel num_threads(num_proc_) 
  //     {
  //       #pragma omp for 
  //       for (time_step=0; time_step<=N_; ++time_step) {
  //         const int robot_id = omp_get_thread_num();
  //         const std::pair<double, double> filter_pair
  //             = split_ocps_[time_step].stageCostAndConstraintsViolation(
  //                 robots_[robot_id], t+time_step*dtau_, dtau_, q_[time_step], 
  //                 v_[time_step], a_[time_step], f_[time_step], u_[time_step]);
  //         costs_.coeffRef(time_step) = filter_pair.first;
  //         constraints_violations_.coeffRef(time_step) = filter_pair.second;
  //         split_ocps_[time_step].getStateDirection(dq_[time_step], 
  //                                                  dv_[time_step]);
  //         if (time_step == N_-1) {
  //           costs_.coeffRef(time_step) 
  //               += split_ocps_[time_step].terminalCost(robots_[robot_id], 
  //                                                      t+T_, q_[N_-1], v_[N_-1]);
  //         }
  //       }
  //     } // #pragma omp parallel num_threads(num_proc_)
  //     filter_.augment(costs_.sum(), constraints_violations_.sum());
  //   }
  //   while (primal_step_size > min_step_size_) {
  //     #pragma omp parallel num_threads(num_proc_) 
  //     {
  //       #pragma omp for 
  //       for (time_step=0; time_step<=N_; ++time_step) {
  //         const int robot_id = omp_get_thread_num();
  //         if (time_step == 0) {
  //           const std::pair<double, double> filter_pair
  //               = split_ocps_[0].stageCostAndConstraintsViolation(
  //                   robots_[robot_id], primal_step_size, t+time_step*dtau_, 
  //                   dtau_, q, v, 
  //                   Eigen::VectorXd::Zero(robots_[robot_id].dimv()), 
  //                   Eigen::VectorXd::Zero(robots_[robot_id].dimv()),
  //                   q_[0], v_[0], a_[0], f_[0], u_[0]);
  //           costs_.coeffRef(0) = filter_pair.first;
  //           constraints_violations_.coeffRef(0) = filter_pair.second;
  //         }
  //         else {
  //           const std::pair<double, double> filter_pair
  //               = split_ocps_[time_step].stageCostAndConstraintsViolation(
  //                   robots_[robot_id], primal_step_size, t+time_step*dtau_, 
  //                   dtau_, q_[time_step-1], v_[time_step-1], dq_[time_step-1], 
  //                   dv_[time_step-1], q_[time_step], v_[time_step], 
  //                   a_[time_step], f_[time_step], u_[time_step]);
  //           costs_.coeffRef(time_step) = filter_pair.first;
  //           constraints_violations_.coeffRef(time_step) = filter_pair.second;
  //         }
  //         if (time_step == N_-1) {
  //           costs_.coeffRef(time_step) 
  //               += split_ocps_[time_step].terminalCost(robots_[robot_id], 
  //                                                      primal_step_size, t+T_, 
  //                                                      q_[N_-1], v_[N_-1]);
  //         }
  //       }
  //     } // #pragma omp parallel num_threads(num_proc_)
  //     const double cost_sum = costs_.sum();
  //     const double constraints_violation_sum = constraints_violations_.sum();
  //     if (filter_.isAccepted(cost_sum, constraints_violation_sum)) {
  //       filter_.augment(cost_sum, constraints_violation_sum);
  //       break;
  //     }
  //     primal_step_size *= step_size_reduction_rate_;
  //   }
  // }  // end if (use_line_search) 
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].updatePrimal(robots_[robot_id], primal_step_size, dtau_, s_[i]);
    split_ocps_[i].updateDual(dual_step_size);
    s_old_[i].lmd = s_[i].lmd;
    s_old_[i].gmm = s_[i].gmm;
    s_old_[i].q = s_[i].q;
    s_old_[i].v = s_[i].v;
    aux_mat_old_[i] = aux_mat_[i];
  }
} 


void ParNMPC::getInitialControlInput(Eigen::VectorXd& u) {
  assert(u.size() == robots_[0].dimv());
  u = s_[0].u;
}


void ParNMPC::getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) {
  // assert(Kq.rows() == robots_[0].dimv());
  // assert(Kq.cols() == robots_[0].dimv());
  // assert(Kv.rows() == robots_[0].dimv());
  // assert(Kv.cols() == robots_[0].dimv());
  // split_ocps_[0].getStateFeedbackGain(Kq, Kv);
}


bool ParNMPC::setStateTrajectory(const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  Eigen::VectorXd q_normalized = q;
  robots_[0].normalizeConfiguration(q_normalized);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    s_[i].v = v;
    s_new_[i].v = v;
    s_old_[i].v = v;
    s_[i].q = q_normalized;
    s_new_[i].q = q_normalized;
    s_old_[i].q = q_normalized;
  }
  bool feasible = isCurrentSolutionFeasible();
  if (feasible) {
    initConstraints();
  }
  return feasible;
}


bool ParNMPC::setStateTrajectory(const Eigen::VectorXd& q0, 
                                 const Eigen::VectorXd& v0, 
                                 const Eigen::VectorXd& qN, 
                                 const Eigen::VectorXd& vN) {
  assert(q0.size() == robots_[0].dimq());
  assert(v0.size() == robots_[0].dimv());
  assert(qN.size() == robots_[0].dimq());
  assert(vN.size() == robots_[0].dimv());
  Eigen::VectorXd q0_normalized = q0;
  robots_[0].normalizeConfiguration(q0_normalized);
  Eigen::VectorXd qN_normalized = qN;
  robots_[0].normalizeConfiguration(qN_normalized);
  const Eigen::VectorXd a = (vN-v0) / N_;
  Eigen::VectorXd dqN = Eigen::VectorXd::Zero(robots_[0].dimv());
  robots_[0].subtractConfiguration(qN, q0, dqN);
  const Eigen::VectorXd v = dqN / N_;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].a = a;
    s_new_[i].a = a;
    s_old_[i].a = a;
    s_[i].v = v0 + i * a;
    s_new_[i].v = s_[i].v;
    s_old_[i].v = s_[i].v;
    s_[i].q = q0;
    robots_[0].integrateConfiguration(v, (double)i, s_[i].q);
    s_new_[i].q = s_[i].q;
    s_old_[i].q = s_[i].q;
  }
  bool feasible = isCurrentSolutionFeasible();
  if (feasible) {
    initConstraints();
  }
  return feasible;
}


void ParNMPC::setContactSequence(
    const std::vector<std::vector<bool>>& contact_sequence) {
  assert(contact_sequence.size() == N_);
  for (int i=0; i<N_; ++i) {
    assert(contact_sequence[i].size() == robots_[0].max_point_contacts());
  }
  contact_sequence_ = contact_sequence;
}


void ParNMPC::resetLineSearchFilter() {
  filter_.clear();
}


double ParNMPC::KKTError(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  Eigen::VectorXd error = Eigen::VectorXd::Zero(N_);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
  const int robot_id = omp_get_thread_num();
  robots_[robot_id].setContactStatus(contact_sequence_[i]);
    if (i == 0) {
      error(i) = split_ocps_[i].squaredKKTErrorNorm(robots_[robot_id], 
                                                    t+i*dtau_, dtau_, q, v, 
                                                    s_[i], s_[i+1].lmd, 
                                                    s_[i+1].gmm, s_[i+1].q,
                                                    s_[i+1].v);
    }
    else if (i<N_-1) {
      error(i) = split_ocps_[i].squaredKKTErrorNorm(robots_[robot_id], 
                                                    t+i*dtau_, dtau_, 
                                                    s_[i-1].q, s_[i-1].v, s_[i], 
                                                    s_[i+1].lmd, s_[i+1].gmm, 
                                                    s_[i+1].q, s_[i+1].v);
    }
  }
  return std::sqrt(error.sum());
}


void ParNMPC::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
    std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
    std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
    std::cout << "f[" << i << "] = " << s_[i].f.transpose() << std::endl;
    std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    std::cout << "mu[" << i << "] = " << s_[i].mu.transpose() << std::endl;
  }
}


bool ParNMPC::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    const bool feasible = split_ocps_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time step " << i << std::endl;
      return false;
    }
  }
  return true;
}


void ParNMPC::initConstraints() {
  for (int i=0; i<N_; ++i) {
    split_ocps_[i].initConstraints(robots_[0], i, dtau_, s_[i]);
  }
}

} // namespace idocp