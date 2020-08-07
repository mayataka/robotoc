#include "idocp/ocp/ocp.hpp"

#include <utility>
#include <cmath>
#include <omp.h>
#include <assert.h>


namespace idocp {

OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunctionInterface>& cost,
         const std::shared_ptr<ConstraintsInterface>& constraints, 
         const double T, const int N, const int num_proc)
  : split_ocps_(N, SplitOCP(robot, cost, constraints)),
    terminal_ocp_(robot, cost, constraints),
    robots_(num_proc, robot),
    filter_(),
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
    f_(N, Eigen::VectorXd::Zero(robot.max_dimf())),
    mu_(N, Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf())),
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
  for (int i=0; i<=N; ++i) {
    robot.normalizeConfiguration(q_[i]);
  }
  bool feasible = isCurrentSolutionFeasible();
  initConstraints();
}


OCP::OCP()
  : split_ocps_(),
    terminal_ocp_(),
    robots_(),
    filter_(),
    T_(0),
    dtau_(0),
    step_size_reduction_rate_(0),
    min_step_size_(0),
    N_(0),
    num_proc_(0),
    q_(),
    v_(),
    a_(),
    u_(),
    beta_(),
    f_(),
    mu_(),
    lmd_(),
    gmm_(),
    dq_(),
    dv_(),
    sq_(),
    sv_(),
    Pqq_(),
    Pqv_(),
    Pvq_(),
    Pvv_(),
    primal_step_sizes_(),
    dual_step_sizes_(),
    costs_(), 
    constraints_violations_(),
    cost_derivative_dot_direction_(),
    contact_sequence_() {
}


OCP::~OCP() {
}


void OCP::solveLQR(const double t, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const bool use_line_search) {
  int time_step;
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for  
    for (time_step=0; time_step<=N_; ++time_step) {
      if (time_step < N_) {
        const int robot_id = omp_get_thread_num();
        robots_[robot_id].setContactStatus(contact_sequence_[time_step]);
        split_ocps_[time_step].linearizeOCP(robots_[robot_id], 
                                            t+time_step*dtau_, dtau_, 
                                            lmd_[time_step], gmm_[time_step],
                                            q_[time_step], v_[time_step], 
                                            a_[time_step], u_[time_step], 
                                            f_[time_step], mu_[time_step],
                                            lmd_[time_step+1], 
                                            gmm_[time_step+1], q_[time_step+1], 
                                            v_[time_step+1]);
      }
      else {
        const int robot_id = omp_get_thread_num();
        terminal_ocp_.linearizeOCP(robots_[robot_id], t+T_, lmd_[N_], gmm_[N_], 
                                   q_[N_], v_[N_], Pqq_[N_], Pqv_[N_], Pvq_[N_], 
                                   Pvv_[N_], sq_[N_], sv_[N_]);
      }
    }
  } // #pragma omp parallel num_threads(num_proc_)
  for (time_step=N_-1; time_step>=0; --time_step) {
    split_ocps_[time_step].backwardRiccatiRecursion(dtau_, Pqq_[time_step+1], 
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
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  robots_[0].subtractConfiguration(q, q_[0], dq_[0]);
  dv_[0] = v - v_[0];
  for (time_step=0; time_step<N_; ++time_step) {
    split_ocps_[time_step].forwardRiccatiRecursion(dtau_, dq_[time_step], 
                                                   dv_[time_step], 
                                                   dq_[time_step+1], 
                                                   dv_[time_step+1]);
  }
  #pragma omp parallel num_threads(num_proc_) 
  {
    #pragma omp for 
    for (time_step=0; time_step<N_; ++time_step) {
      split_ocps_[time_step].computeCondensedDirection(dtau_, dq_[time_step], 
                                                       dv_[time_step]);
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
          if (time_step < N_) {
            const int robot_id = omp_get_thread_num();
            const std::pair<double, double> filter_pair
                = split_ocps_[time_step].costAndConstraintsViolation(
                    robots_[robot_id], t+time_step*dtau_, dtau_, q_[time_step], 
                    v_[time_step], a_[time_step], u_[time_step], f_[time_step]);
            costs_.coeffRef(time_step) = filter_pair.first;
            constraints_violations_.coeffRef(time_step) = filter_pair.second;
          }
          else {
            const int robot_id = omp_get_thread_num();
            costs_.coeffRef(N_) = terminal_ocp_.terminalCost(robots_[robot_id], 
                                                             t+T_, q_[N_], 
                                                             v_[N_]);
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
                = split_ocps_[time_step].costAndConstraintsViolation(
                    robots_[robot_id], primal_step_size, t+time_step*dtau_, 
                    dtau_, q_[time_step], v_[time_step], a_[time_step], 
                    u_[time_step], f_[time_step], q_[time_step+1], 
                    v_[time_step+1], dq_[time_step], dv_[time_step], 
                    dq_[time_step+1], dv_[time_step+1]);
            costs_.coeffRef(time_step) = filter_pair.first;
            constraints_violations_.coeffRef(time_step) = filter_pair.second;
          }
          else {
            const int robot_id = omp_get_thread_num();
            costs_.coeffRef(N_) = terminal_ocp_.terminalCost(robots_[robot_id], 
                                                             primal_step_size, 
                                                             t+T_, q_[N_], v_[N_], 
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
        split_ocps_[time_step].updatePrimal(robots_[robot_id], primal_step_size, 
                                            dtau_, Pqq_[time_step], 
                                            Pqv_[time_step], Pvq_[time_step], 
                                            Pvv_[time_step], sq_[time_step], 
                                            sv_[time_step], dq_[time_step], 
                                            dv_[time_step], lmd_[time_step], 
                                            gmm_[time_step], q_[time_step], 
                                            v_[time_step], a_[time_step], 
                                            u_[time_step], beta_[time_step], 
                                            f_[time_step], mu_[time_step]);
        split_ocps_[time_step].updateDual(dual_step_size);
      }
      else {
        const int robot_id = omp_get_thread_num();
        terminal_ocp_.updatePrimal(robots_[robot_id], primal_step_size, 
                                   Pqq_[N_], Pqv_[N_], Pvq_[N_], Pvv_[N_], 
                                   sq_[N_], sv_[N_], dq_[N_], dv_[N_], 
                                   lmd_[N_], gmm_[N_], q_[N_], v_[N_]);
        terminal_ocp_.updateDual(dual_step_size);
      }
    }
  } // #pragma omp parallel num_threads(num_proc_)
} 


void OCP::getInitialControlInput(Eigen::VectorXd& u) {
  assert(u.size() == robots_[0].dimv());
  u = u_[0];
}


void OCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) {
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  split_ocps_[0].getStateFeedbackGain(Kq, Kv);
}


bool OCP::setStateTrajectory(const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  for (int i=0; i<=N_; ++i) {
    v_[i] = v;
  }
  for (int i=0; i<=N_; ++i) {
    q_[i] = q;
    robots_[0].normalizeConfiguration(q_[i]);
  }
  bool feasible = isCurrentSolutionFeasible();
  if (feasible) {
    initConstraints();
  }
  return feasible;
}


bool OCP::setStateTrajectory(const Eigen::VectorXd& q0, 
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
  }
  for (int i=0; i<=N_; ++i) {
    q_[i] = q0;
    robots_[0].integrateConfiguration(v, (double)i, q_[i]);
    robots_[0].normalizeConfiguration(q_[i]);
  }
  bool feasible = isCurrentSolutionFeasible();
  if (feasible) {
    initConstraints();
  }
  return feasible;
}


void OCP::setContactSequence(
    const std::vector<std::vector<bool>>& contact_sequence) {
  assert(contact_sequence.size() == N_);
  for (int i=0; i<N_; ++i) {
    assert(contact_sequence[i].size() == robots_[0].max_point_contacts());
  }
  contact_sequence_ = contact_sequence;
}


void OCP::resetLineSearchFilter() {
  filter_.clear();
}


double OCP::KKTError(const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v) {
  assert(q.size() == q_[0].size());
  assert(v.size() == v_[0].size());
  double error = 0;
  error += (q-q_[0]).squaredNorm();
  error += (v-v_[0]).squaredNorm();
  for (int i=0; i<N_; ++i) {
    error += split_ocps_[i].squaredKKTErrorNorm(robots_[0], t+i*dtau_, dtau_, 
                                                lmd_[i], gmm_[i], q_[i], v_[i], 
                                                a_[i], u_[i], beta_[i], f_[i], 
                                                mu_[i], lmd_[i+1], gmm_[i+1], 
                                                q_[i+1], v_[i+1]);
  }
  error += terminal_ocp_.squaredKKTErrorNorm(robots_[0], t+T_, 
                                             lmd_[N_], gmm_[N_], q_[N_], v_[N_]);
  return std::sqrt(error);
}


void OCP::printSolution() const {
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


bool OCP::isCurrentSolutionFeasible() {
  for (int time_step=0; time_step<N_; ++time_step) {
    const bool feasible = split_ocps_[time_step].isFeasible(robots_[0], 
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
    split_ocps_[time_step].initConstraints(robots_[0], time_step, dtau_, 
                                           q_[time_step], v_[time_step], 
                                           a_[time_step], u_[time_step]);
  }
}

} // namespace idocp