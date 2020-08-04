#include "idocp/ocp/ocp.hpp"

#include <cmath>
#include <assert.h>
#include <omp.h>


namespace idocp {

OCP::OCP(const Robot& robot, const CostFunctionInterface& cost,
         const ConstraintsInterface& constraints, const double T, const int N, 
         const int num_proc)
  : split_OCPs_(N, SplitOCP(robot, cost, constraints)),
    split_terminal_OCP_(robot, cost, constraints),
    robots_(num_proc, robot),
    filter_(),
    cost_(cost),
    constraints_(cost),
    T_(T),
    dtau_(T/N),
    step_size_reduction_rate_(0.75),
    min_step_size_(0.05),
    N_(N),
    num_proc_(num_proc),
    q_(N+1, Eigen::VectorXd::Zero(robot.dimq())),
    v_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    a_(N, Eigen::VectorXd::Zero(robot.dimv())),
    u_(N, Eigen::VectorXd::Zero(robot.dimv())),
    beta_(N, Eigen::VectorXd::Zero(robot.dimv())),
    lmd_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    gmm_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    dq_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    dv_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    sq_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    sv_(N+1, Eigen::VectorXd::Zero(robot.dimv())),
    Pqq_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pqv_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pvq_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Pvv_(N+1, Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    primal_step_sizes_(Eigen::VectorXd::Zero(N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N)),
    costs_(Eigen::VectorXd::Zero(N+1)), 
    constraints_violations_(Eigen::VectorXd::Zero(N)),
    cost_derivative_dot_direction_(Eigen::VectorXd::Zero(N+1)),
    contact_sequence_(N, std::vector<bool>(robot.max_point_contacts(), false)) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  bool feasible = isCurrentSolutionFeasible();
  initConstraints();
  activateAllContacts();
}


OCP::~OCP() {
}


void OCP::solveSQP(const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const bool use_line_search) {
  int time_step;
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for  
    for (time_step=0; time_step<=N_; ++time_step) {
      if (time_step < N_) {
        const int robot_id = omp_get_thread_num();
        robots_[robot_id].setActiveContacts(contact_sequence_[time_step]);
        split_OCPs_[time_step].linearizeOCP(robots_[robot_id], 
                                            t+time_step*dtau_, dtau_, 
                                            lmd_[time_step], gmm_[time_step],
                                            q_[time_step], v_[time_step], 
                                            a_[time_step], u_[time_step], 
                                            lmd_[time_step+1], 
                                            gmm_[time_step+1], q_[time_step+1], 
                                            v_[time_step+1]);
      }
      else {
        const int robot_id = omp_get_thread_num();
        split_terminal_OCP_.linearizeOCP(robots_[robot_id], t+T_, lmd_[N_], 
                                         gmm_[N_], q_[N_], v_[N_], Pqq_[N_], 
                                         Pqv_[N_], Pvq_[N_], Pvv_[N_], sq_[N_], 
                                         sv_[N_]);
      }
    }
  } // #pragma omp parallel num_threads(num_proc_)
  for (time_step=N_-1; time_step>=0; --time_step) {
    split_OCPs_[time_step].backwardRiccatiRecursion(dtau_, Pqq_[time_step+1], 
                                                    Pqv_[time_step+1], 
                                                    Pvq_[time_step+1], 
                                                    Pvv_[time_step+1], 
                                                    sq_[time_step+1], 
                                                    sv_[time_step+1], 
                                                    Pqq_[time_step], 
                                                    Pqv_[time_step], 
                                                    Pvq_[time_step], 
                                                    Pvv_[time_step], 
                                                    sq_[time_step], 
                                                    sv_[time_step]);
  }
  assert(q.size() == q_[0].size());
  assert(v.size() == v_[0].size());
  dq_[0] = q - q_[0];
  dv_[0] = v - v_[0];
  for (time_step=0; time_step<N_; ++time_step) {
    split_OCPs_[time_step].forwardRiccatiRecursion(dtau_, dq_[time_step], 
                                                   dv_[time_step], 
                                                   dq_[time_step+1], 
                                                   dv_[time_step+1]);
  }
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for 
    for (time_step=0; time_step<N_; ++time_step) {
      const int robot_id = omp_get_thread_num();
      split_OCPs_[time_step].computeCondensedDirection(robots_[robot_id], 
                                                       dtau_, dq_[time_step], 
                                                       dv_[time_step]);
      primal_step_sizes_.coeffRef(time_step) 
          = split_OCPs_[time_step].maxPrimalStepSize();
      dual_step_sizes_.coeffRef(time_step) 
          = split_OCPs_[time_step].maxDualStepSize();
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
          if (time_step < N_) {
            const int robot_id = omp_get_thread_num();
            const std::pair<double, double> filter_pair
                = split_OCPs_[time_step].costAndConstraintsViolation(
                    robots_[robot_id], t+time_step*dtau_, dtau_, q_[time_step], 
                    v_[time_step], a_[time_step], u_[time_step]);
            costs_.coeffRef(time_step) = filter_pair.first;
            constraints_violations_.coeffRef(time_step) = filter_pair.second;
          }
          else {
            const int robot_id = omp_get_thread_num();
            costs_.coeffRef(N_) = split_terminal_OCP_.terminalCost(
                robots_[robot_id], t+T_, q_[N_], v_[N_]);
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
          if (time_step < N_) {
            const int robot_id = omp_get_thread_num();
            const std::pair<double, double> filter_pair
                = split_OCPs_[time_step].costAndConstraintsViolation(
                    robots_[robot_id], primal_step_size, t+time_step*dtau_, 
                    dtau_, q_[time_step], v_[time_step], a_[time_step], 
                    u_[time_step], q_[time_step+1], v_[time_step+1], 
                    dq_[time_step], dv_[time_step], dq_[time_step+1], 
                    dv_[time_step+1]);
            costs_.coeffRef(time_step) = filter_pair.first;
            constraints_violations_.coeffRef(time_step) = filter_pair.second;
          }
          else {
            const int robot_id = omp_get_thread_num();
            costs_.coeffRef(N_) = split_OCPs_[N_].terminalCost(
                robots_[robot_id], primal_step_size, t+T_, q_[N_], v_[N_], 
                dq_[N_], dv_[N_]);
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
    for (time_step=0; time_step<=N_; ++time_step) {
      if (time_step < N_) {
        const int robot_id = omp_get_thread_num();
        split_OCPs_[time_step].updateDual(dual_step_size);
        split_OCPs_[time_step].updatePrimal(robots_[robot_id], primal_step_size, 
                                            dtau_, dq_[time_step], 
                                            dv_[time_step], Pqq_[time_step], 
                                            Pqv_[time_step], Pvq_[time_step], 
                                            Pvv_[time_step], sq_[time_step], 
                                            sv_[time_step], q_[time_step], 
                                            v_[time_step], a_[time_step], 
                                            u_[time_step], beta_[time_step], 
                                            lmd_[time_step], gmm_[time_step]);
      }
      else {
        const int robot_id = omp_get_thread_num();
        split_terminal_OCP_.updatePrimal(robots_[robot_id], primal_step_size, 
                                         dq_[N_], dv_[N_], Pqq_[N_], Pqv_[N_], 
                                         Pvq_[N_], Pvv_[N_], sq_[N_], sv_[N_], 
                                         q_[N_], v_[N_], lmd_[N_], gmm_[N_]);
      }
    }
  } // #pragma omp parallel num_threads(num_proc_)
} 


void OCP::getInitialControlInput(Eigen::VectorXd& u) {
  assert(u.size() == u_[0].size());
  u = u_[0];
}


void OCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) {
  assert(Kq.cols() == Kq.rows());
  assert(Kv.cols() == Kv.rows());
  split_OCPs_[0].getStateFeedbackGain(Kq, Kv);
}


void OCP::setStateTrajectory(const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(q.size() == q_[0].size());
  assert(v.size() == v_[0].size());
  for (int i=0; i<=N_; ++i) {
    v_[i] = v;
  }
  for (int i=0; i<=N_; ++i) {
    q_[i] = q;
  }
  bool feasible = isCurrentSolutionFeasible();
  initConstraints();
}


void OCP::setStateTrajectory(const Eigen::VectorXd& q0, 
                             const Eigen::VectorXd& v0, 
                             const Eigen::VectorXd& qN, 
                             const Eigen::VectorXd& vN) {
  assert(q0.size() == q_[0].size());
  assert(v0.size() == v_[0].size());
  assert(qN.size() == q_[0].size());
  assert(vN.size() == v_[0].size());
  const Eigen::VectorXd a = (vN-v0) / N_;
  const Eigen::VectorXd v = (qN-q0) / N_;
  for (int i=0; i<N_; ++i) {
    a_[i] = a;
  }
  for (int i=0; i<=N_; ++i) {
    v_[i] = v0 + i * a;
  }
  for (int i=0; i<=N_; ++i) {
    q_[i] = q0 + i * v;
  }
  bool feasible = isCurrentSolutionFeasible();
  initConstraints();
}


double OCP::KKTError(const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v) {
  assert(q.size() == q_[0].size());
  assert(v.size() == v_[0].size());
  double error = 0;
  error += (q-q_[0]).squaredNorm();
  error += (v-v_[0]).squaredNorm();
  for (int i=0; i<N_; ++i) {
    error += split_OCPs_[i].squaredKKTErrorNorm(robots_[0], t+i*dtau_, dtau_, 
                                                lmd_[i], gmm_[i], q_[i], v_[i], 
                                                a_[i], u_[i], beta_[i], 
                                                lmd_[i+1], gmm_[i+1], 
                                                q_[i+1], v_[i+1]);
  }
  error += split_terminal_OCP_.squaredKKTErrorNorm(robots_[0], t+T_, lmd_[N_], 
                                                   gmm_[N_], q_[N_], v_[N_]);
  return std::sqrt(error);
}


void OCP::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "q[" << i << "] = " << q_[i].transpose() << std::endl;
    std::cout << "v[" << i << "] = " << v_[i].transpose() << std::endl;
    std::cout << "a[" << i << "] = " << a_[i].transpose() << std::endl;
    std::cout << "u[" << i << "] = " << u_[i].transpose() << std::endl;
  }
  std::cout << "q[" << N_ << "] = " << q_[N_].transpose() << std::endl;
  std::cout << "v[" << N_ << "] = " << v_[N_].transpose() << std::endl;
}


bool OCP::isCurrentSolutionFeasible() {
  for (int time_step=0; time_step<N_; ++time_step) {
    const bool feasible = split_OCPs_[time_step].isFeasible(robots_[0], 
                                                            q_[time_step], 
                                                            v_[time_step], 
                                                            a_[time_step], 
                                                            u_[time_step]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time step " << time_step << std::endl;
      return false;
    }
  }
  return true;
}


void OCP::initConstraints() {
  for (int time_step=0; time_step<N_; ++time_step) {
    split_OCPs_[time_step].initConstraints(robots_[0], time_step, dtau_, 
                                           q_[time_step], v_[time_step], 
                                           a_[time_step], u_[time_step]);
  }
}

void OCP::activateAllContacts() {
  for (int i=0; i<contact_sequence_.size(); ++i) {
    for (int j=0; j<robots_[0].max_point_contacts(); ++j) {
      contact_sequence_[i][j] = true;
    }
  }
}

} // namespace idocp