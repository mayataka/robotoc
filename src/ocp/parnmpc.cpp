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
    lmd_(N, Eigen::VectorXd::Zero(robot.dimv())),
    gmm_(N, Eigen::VectorXd::Zero(robot.dimv())),
    mu_(N, Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf())),
    a_(N, Eigen::VectorXd::Zero(robot.dimv())),
    f_(N, Eigen::VectorXd::Zero(robot.max_dimf())),
    q_(N+1, Eigen::VectorXd::Zero(robot.dimq())),
    v_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    u_(N, Eigen::VectorXd::Zero(robot.dimv())),
    beta_(N, Eigen::VectorXd::Zero(robot.dimv())),
    lmd_old_(N, Eigen::VectorXd::Zero(robot.dimv())),
    gmm_old_(N, Eigen::VectorXd::Zero(robot.dimv())),
    q_old_(N, Eigen::VectorXd::Zero(robot.dimq())),
    v_old_(N, Eigen::VectorXd::Zero(robot.dimv())),
    dq_(N, Eigen::VectorXd::Zero(robot.dimv())),
    dv_(N, Eigen::VectorXd::Zero(robot.dimv())),
    aux_mat_(N, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    aux_mat_old_(N+1, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    primal_step_sizes_(Eigen::VectorXd::Zero(N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N)),
    costs_(Eigen::VectorXd::Zero(N)), 
    constraints_violations_(Eigen::VectorXd::Zero(N)),
    contact_sequence_(N, std::vector<bool>(robot.max_point_contacts(), false)) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  for (int i=0; i<=N; ++i) {
    robot.normalizeConfiguration(q_[i]);
    robot.normalizeConfiguration(q_old_[i]);
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
    lmd_(),
    gmm_(),
    mu_(),
    a_(),
    f_(),
    q_(),
    v_(),
    u_(),
    beta_(),
    lmd_old_(),
    gmm_old_(),
    q_old_(),
    v_old_(),
    dq_(),
    dv_(),
    aux_mat_(),
    aux_mat_old_(),
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
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for  
    for (time_step=0; time_step<N_; ++time_step) {
      const int robot_id = omp_get_thread_num();
      robots_[robot_id].setContactStatus(contact_sequence_[time_step]);
      bool is_terminal_ocp;
      if (time_step == 0) {
        split_ocps_[time_step].coarseUpdate(robots_[robot_id], 
                                            t+time_step*dtau_, dtau_, q, v, 
                                            lmd_[time_step], gmm_[time_step],
                                            mu_[time_step], a_[time_step], 
                                            f_[time_step], q_[time_step], 
                                            v_[time_step], u_[time_step],
                                            lmd_[time_step+1], 
                                            gmm_[time_step+1], q_[time_step+1], 
                                            v_[time_step+1], 
                                            aux_mat_old_[time_step+1], 
                                            aux_mat_[time_step], false);
      }
      else if (time_step < N_-1) {
        split_ocps_[time_step].coarseUpdate(robots_[robot_id], 
                                            t+time_step*dtau_, dtau_, 
                                            q_[time_step-1], v_[time_step-1], 
                                            lmd_[time_step], gmm_[time_step],
                                            mu_[time_step], a_[time_step], 
                                            f_[time_step], q_[time_step], 
                                            v_[time_step], u_[time_step],
                                            lmd_[time_step+1], 
                                            gmm_[time_step+1], q_[time_step+1], 
                                            v_[time_step+1], 
                                            aux_mat_old_[time_step+1], 
                                            aux_mat_[time_step], false);
      }
      else {
        split_ocps_[time_step].computeTerminalCostDerivatives(
            robots_[robot_id], t+T_, q_[N_-1], v_[N_-1], lmd_[N_], gmm_[N_]);
        split_ocps_[time_step].coarseUpdate(robots_[robot_id], 
                                            t+time_step*dtau_, dtau_, 
                                            q_[time_step-1], v_[time_step-1], 
                                            lmd_[time_step], gmm_[time_step],
                                            mu_[time_step], a_[time_step], 
                                            f_[time_step], q_[time_step], 
                                            v_[time_step], u_[time_step],
                                            lmd_[time_step+1], 
                                            gmm_[time_step+1], q_[time_step+1], 
                                            v_[time_step+1], 
                                            aux_mat_old_[time_step+1], 
                                            aux_mat_[time_step], true);
      }
    }
  } // #pragma omp parallel num_threads(num_proc_)
  for (time_step=N_-1; time_step>=0; --time_step) {
    split_ocps_[time_step].backwardCollectionSerial(lmd_old_[time_step+1], 
                                                    gmm_old_[time_step+1],
                                                    lmd_[time_step+1], 
                                                    gmm_[time_step+1],
                                                    lmd_[time_step], 
                                                    gmm_[time_step]);
  }
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for  
    for (time_step=0; time_step<N_; ++time_step) {
      const int robot_id = omp_get_thread_num();
      robots_[robot_id].setContactStatus(contact_sequence_[time_step]);
      split_ocps_[time_step].backwardCollectionParallel(robots_[robot_id]);
    }
  } // #pragma omp parallel num_threads(num_proc_)
  split_ocps_[0].forwardCollectionSerial(robots_[0], q, v, q, v, q_[0], v_[0]);
  for (time_step=1; time_step<N_; ++time_step) {
    split_ocps_[time_step].forwardCollectionSerial(robots_[0], 
                                                   q_old_[time_step-1], 
                                                   v_old_[time_step-1], 
                                                   q_[time_step-1], 
                                                   v_[time_step-1], 
                                                   q_[time_step], 
                                                   v_[time_step]);
  }
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for 
    for (time_step=0; time_step<N_; ++time_step) {
      const int robot_id = omp_get_thread_num();
      robots_[robot_id].setContactStatus(contact_sequence_[time_step]);
      split_ocps_[time_step].forwardCollectionParallel(robots_[robot_id]);
      split_ocps_[time_step].computePrimalAndDualDirection(robots_[robot_id], 
                                                           dtau_, 
                                                           lmd_[time_step],
                                                           gmm_[time_step],
                                                           mu_[time_step],
                                                           a_[time_step],
                                                           f_[time_step],
                                                           q_[time_step],
                                                           v_[time_step],
                                                           u_[time_step]);
      primal_step_sizes_.coeffRef(time_step) 
          = split_ocps_[time_step].maxPrimalStepSize();
      dual_step_sizes_.coeffRef(time_step) 
          = split_ocps_[time_step].maxDualStepSize();
    }
  } // #pragma omp parallel num_threads(num_proc_)
  double primal_step_size = primal_step_sizes_.minCoeff();
  const double dual_step_size = dual_step_sizes_.minCoeff();
  if (use_line_search) {
    // If filter is empty, augment the current solution to the filter.
    if (filter_.isEmpty()) {
      #pragma omp parallel num_threads(num_proc_) 
      {
        #pragma omp for 
        for (time_step=0; time_step<=N_; ++time_step) {
          const int robot_id = omp_get_thread_num();
          const std::pair<double, double> filter_pair
              = split_ocps_[time_step].stageCostAndConstraintsViolation(
                  robots_[robot_id], t+time_step*dtau_, dtau_, q_[time_step], 
                  v_[time_step], a_[time_step], f_[time_step], u_[time_step]);
          costs_.coeffRef(time_step) = filter_pair.first;
          constraints_violations_.coeffRef(time_step) = filter_pair.second;
          split_ocps_[time_step].getStateDirection(dq_[time_step], 
                                                   dv_[time_step]);
          if (time_step == N_-1) {
            costs_.coeffRef(time_step) 
                += split_ocps_[time_step].terminalCost(robots_[robot_id], 
                                                       t+T_, q_[N_-1], v_[N_-1]);
          }
        }
      } // #pragma omp parallel num_threads(num_proc_)
      filter_.augment(costs_.sum(), constraints_violations_.sum());
    }
    while (primal_step_size > min_step_size_) {
      #pragma omp parallel num_threads(num_proc_) 
      {
        #pragma omp for 
        for (time_step=0; time_step<=N_; ++time_step) {
          const int robot_id = omp_get_thread_num();
          if (time_step == 0) {
            const std::pair<double, double> filter_pair
                = split_ocps_[0].stageCostAndConstraintsViolation(
                    robots_[robot_id], primal_step_size, t+time_step*dtau_, 
                    dtau_, q, v, 
                    Eigen::VectorXd::Zero(robots_[robot_id].dimv()), 
                    Eigen::VectorXd::Zero(robots_[robot_id].dimv()),
                    q_[0], v_[0], a_[0], f_[0], u_[0]);
            costs_.coeffRef(0) = filter_pair.first;
            constraints_violations_.coeffRef(0) = filter_pair.second;
          }
          else {
            const std::pair<double, double> filter_pair
                = split_ocps_[time_step].stageCostAndConstraintsViolation(
                    robots_[robot_id], primal_step_size, t+time_step*dtau_, 
                    dtau_, q_[time_step-1], v_[time_step-1], dq_[time_step-1], 
                    dv_[time_step-1], q_[time_step], v_[time_step], 
                    a_[time_step], f_[time_step], u_[time_step]);
            costs_.coeffRef(time_step) = filter_pair.first;
            constraints_violations_.coeffRef(time_step) = filter_pair.second;
          }
          if (time_step == N_-1) {
            costs_.coeffRef(time_step) 
                += split_ocps_[time_step].terminalCost(robots_[robot_id], 
                                                       primal_step_size, t+T_, 
                                                       q_[N_-1], v_[N_-1]);
          }
        }
      } // #pragma omp parallel num_threads(num_proc_)
      const double cost_sum = costs_.sum();
      const double constraints_violation_sum = constraints_violations_.sum();
      if (filter_.isAccepted(cost_sum, constraints_violation_sum)) {
        filter_.augment(cost_sum, constraints_violation_sum);
        break;
      }
      primal_step_size *= step_size_reduction_rate_;
    }
  }  // end if (use_line_search) 
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for 
    for (time_step=0; time_step<N_; ++time_step) {
      const int robot_id = omp_get_thread_num();
      split_ocps_[time_step].updatePrimal(robots_[robot_id], primal_step_size, 
                                          dtau_, lmd_[time_step], 
                                          gmm_[time_step], a_[time_step], 
                                          f_[time_step], mu_[time_step], 
                                          q_[time_step], v_[time_step], 
                                          u_[time_step], beta_[time_step]);
      split_ocps_[time_step].updateDual(dual_step_size);
      lmd_old_[time_step] = lmd_[time_step];
      gmm_old_[time_step] = gmm_[time_step];
      q_old_[time_step] = q_[time_step];
      v_old_[time_step] = v_[time_step];
      aux_mat_old_[time_step] = aux_mat_[time_step];
    }
  } // #pragma omp parallel num_threads(num_proc_)
} 


void ParNMPC::getInitialControlInput(Eigen::VectorXd& u) {
  assert(u.size() == robots_[0].dimv());
  u = u_[0];
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
  for (int i=0; i<=N_; ++i) {
    v_[i] = v;
    v_old_[i] = v;
  }
  for (int i=0; i<=N_; ++i) {
    q_[i] = q;
    robots_[0].normalizeConfiguration(q_[i]);
    q_old_[i] = q_[i];
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
  const Eigen::VectorXd a = (vN-v0) / N_;
  Eigen::VectorXd dqN = Eigen::VectorXd::Zero(robots_[0].dimv());
  robots_[0].subtractConfiguration(qN, q0, dqN);
  const Eigen::VectorXd v = dqN / N_;
  for (int i=0; i<N_; ++i) {
    a_[i] = a;
  }
  for (int i=0; i<=N_; ++i) {
    v_[i] = v0 + i * a;
    v_old_[i] = v0 + i * a;
  }
  for (int i=0; i<=N_; ++i) {
    q_[i] = q0;
    robots_[0].integrateConfiguration(v, (double)i, q_[i]);
    robots_[0].normalizeConfiguration(q_[i]);
    q_old_[i] = q_[i];
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
  assert(q.size() == q_[0].size());
  assert(v.size() == v_[0].size());
  double error = 0;
  error += (q-q_[0]).squaredNorm();
  error += (v-v_[0]).squaredNorm();
  for (int i=0; i<N_; ++i) {
    if (i == N_-1) {
      split_ocps_[i].computeTerminalCostDerivatives(robots_[0], t+T_, q_[i], 
                                                    v_[i], lmd_[i+1], gmm_[i+1]);
    }
    error += split_ocps_[i].squaredKKTErrorNorm(robots_[0], t+i*dtau_, dtau_, 
                                                lmd_[i], gmm_[i], q_[i], v_[i], 
                                                a_[i], u_[i], beta_[i], f_[i], 
                                                mu_[i], lmd_[i+1], gmm_[i+1], 
                                                q_[i+1], v_[i+1]);
  }
  return std::sqrt(error);
}


void ParNMPC::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "q[" << i << "] = " << q_[i].transpose() << std::endl;
    std::cout << "v[" << i << "] = " << v_[i].transpose() << std::endl;
    std::cout << "a[" << i << "] = " << a_[i].transpose() << std::endl;
    std::cout << "u[" << i << "] = " << u_[i].transpose() << std::endl;
    std::cout << "f[" << i << "] = " << f_[i].transpose() << std::endl;
    std::cout << "mu[" << i << "] = " << mu_[i].transpose() << std::endl;
  }
  std::cout << "q[" << N_ << "] = " << q_[N_].transpose() << std::endl;
  std::cout << "v[" << N_ << "] = " << v_[N_].transpose() << std::endl;
}


bool ParNMPC::isCurrentSolutionFeasible() {
  for (int time_step=0; time_step<N_; ++time_step) {
    const bool feasible = split_ocps_[time_step].isFeasible(robots_[0], 
                                                            q_[time_step], 
                                                            v_[time_step], 
                                                            a_[time_step], 
                                                            f_[time_step], 
                                                            u_[time_step]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time step " << time_step << std::endl;
      return false;
    }
  }
  return true;
}


void ParNMPC::initConstraints() {
  for (int time_step=0; time_step<N_; ++time_step) {
    split_ocps_[time_step].initConstraints(robots_[0], time_step, dtau_, 
                                           q_[time_step], v_[time_step], 
                                           a_[time_step], f_[time_step], 
                                           u_[time_step]);
  }
}

} // namespace idocp