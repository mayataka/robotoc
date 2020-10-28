#include "idocp/ocp/ocp.hpp"

#include <utility>
#include <cmath>
#include <omp.h>
#include <assert.h>


namespace idocp {

OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
         const std::shared_ptr<Constraints>& constraints, const double T, 
         const int N, const int num_proc)
  : split_ocps_(N, SplitOCP(robot, cost, constraints)),
    terminal_ocp_(robot, cost, constraints),
    robots_(num_proc, robot),
    contact_sequence_(robot, N),
    filter_(),
    T_(T),
    dtau_(T/N),
    step_size_reduction_rate_(0.75),
    min_step_size_(0.05),
    N_(N),
    num_proc_(num_proc),
    s_(N+1, SplitSolution(robot)),
    d_(N+1, SplitDirection(robot)),
    riccati_(N+1, RiccatiFactorization(robot)),
    primal_step_sizes_(Eigen::VectorXd::Zero(N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N)),
    costs_(Eigen::VectorXd::Zero(N+1)), 
    violations_(Eigen::VectorXd::Zero(N)) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N; ++i) {
    robot.normalizeConfiguration(s_[i].q);
  }
  initConstraints();
}


OCP::OCP() 
  : split_ocps_(),
    terminal_ocp_(),
    robots_(),
    contact_sequence_(),
    filter_(),
    T_(),
    dtau_(),
    step_size_reduction_rate_(),
    min_step_size_(),
    N_(),
    num_proc_(),
    s_(),
    d_(),
    riccati_(),
    primal_step_sizes_(),
    dual_step_sizes_(),
    costs_(), 
    violations_() {
}


OCP::~OCP() {
}


void OCP::updateSolution(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v, const bool use_line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    if (i == 0) {
      const int robot_id = omp_get_thread_num();
      split_ocps_[i].linearizeOCP(robots_[robot_id], 
                                  contact_sequence_.contactStatus(i), t+i*dtau_, 
                                  dtau_, q, s_[i], s_[i+1]);
    }
    else if (i < N_) {
      const int robot_id = omp_get_thread_num();
      split_ocps_[i].linearizeOCP(robots_[robot_id], 
                                  contact_sequence_.contactStatus(i), t+i*dtau_, 
                                  dtau_, s_[i-1].q, s_[i], s_[i+1]);
    }
    else {
      const int robot_id = omp_get_thread_num();
      terminal_ocp_.linearizeOCP(robots_[robot_id], t+T_, s_[N_], riccati_[N_]);
    }
  }
  for (int i=N_-1; i>=0; --i) {
    split_ocps_[i].backwardRiccatiRecursion(dtau_, riccati_[i+1], riccati_[i]);
  }
  robots_[0].subtractConfiguration(q, s_[0].q, d_[0].dq());
  d_[0].dv() = v - s_[0].v;
  for (int i=0; i<N_; ++i) {
    split_ocps_[i].forwardRiccatiRecursion(dtau_, d_[i], d_[i+1]);
  }
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].computeCondensedDirection(robots_[robot_id], dtau_, 
                                             s_[i], d_[i+1], d_[i]);
    primal_step_sizes_.coeffRef(i) = split_ocps_[i].maxPrimalStepSize();
    dual_step_sizes_.coeffRef(i) = split_ocps_[i].maxDualStepSize();
  }
  double primal_step_size = primal_step_sizes_.minCoeff();
  const double dual_step_size = dual_step_sizes_.minCoeff();
  if (use_line_search) {
    // If filter is empty, augment the current solution to the filter.
    if (filter_.isEmpty()) {
      #pragma omp parallel for num_threads(num_proc_)
      for (int i=0; i<=N_; ++i) {
        if (i < N_) {
          const int robot_id = omp_get_thread_num();
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndConstraintViolation(robots_[robot_id], 
                                                          t+i*dtau_, dtau_, 
                                                          s_[i]);
          costs_.coeffRef(i) = filter_pair.first;
          violations_.coeffRef(i) = filter_pair.second;
        }
        else {
          const int robot_id = omp_get_thread_num();
          costs_.coeffRef(N_) = terminal_ocp_.terminalCost(robots_[robot_id], 
                                                           t+T_, s_[N_]);
        }
      }
      filter_.augment(costs_.sum(), violations_.sum());
    }
    while (primal_step_size > min_step_size_) {
      #pragma omp parallel for num_threads(num_proc_)
      for (int i=0; i<=N_; ++i) {
        if (i < N_) {
          const int robot_id = omp_get_thread_num();
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndConstraintViolation(
                  robots_[robot_id], contact_sequence_.contactStatus(i), 
                  primal_step_size, t+i*dtau_, dtau_, s_[i], d_[i], 
                  s_[i+1], d_[i+1]);
          costs_.coeffRef(i) = filter_pair.first;
          violations_.coeffRef(i) = filter_pair.second;
        }
        else {
          const int robot_id = omp_get_thread_num();
          costs_.coeffRef(N_) = terminal_ocp_.terminalCost(robots_[robot_id], 
                                                           primal_step_size, 
                                                           t+T_, s_[N_], d_[N_]);
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
  }  // end if (use_line_search) 
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      const int robot_id = omp_get_thread_num();
      split_ocps_[i].updatePrimal(robots_[robot_id], primal_step_size, dtau_, 
                                  riccati_[i], d_[i], s_[i]);
      split_ocps_[i].updateDual(dual_step_size);
    }
    else {
      const int robot_id = omp_get_thread_num();
      terminal_ocp_.updatePrimal(robots_[robot_id], primal_step_size, 
                                 riccati_[N_], d_[N_], s_[N_]);
      terminal_ocp_.updateDual(dual_step_size);
    }
  }
} 


void OCP::getControlInput(const int stage, Eigen::VectorXd& u) const {
  assert(stage >= 0);
  assert(stage < N_);
  assert(u.size() == robots_[0].dimv());
  u = s_[stage].u;
}


void OCP::getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                               Eigen::MatrixXd& Kv) const {
  assert(stage >= 0);
  assert(stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  split_ocps_[stage].getStateFeedbackGain(Kq, Kv);
}


bool OCP::setStateTrajectory(const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  Eigen::VectorXd q_normalized = q;
  robots_[0].normalizeConfiguration(q_normalized);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    s_[i].v = v;
    s_[i].q = q_normalized;
  }
  initConstraints();
  const bool feasible = isCurrentSolutionFeasible();
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
  Eigen::VectorXd q0_normalized = q0;
  robots_[0].normalizeConfiguration(q0_normalized);
  Eigen::VectorXd qN_normalized = qN;
  robots_[0].normalizeConfiguration(qN_normalized);
  const Eigen::VectorXd a = (vN-v0) / N_;
  Eigen::VectorXd dqN = Eigen::VectorXd::Zero(robots_[0].dimv());
  robots_[0].subtractConfiguration(qN, q0, dqN);
  const Eigen::VectorXd v = dqN / N_;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    s_[i].a = a;
    s_[i].v = v0 + i * a;
    robots_[0].integrateConfiguration(q0, v, (double)i, s_[i].q);
  }
  initConstraints();
  const bool feasible = isCurrentSolutionFeasible();
  return feasible;
}


void OCP::activateContact(const int contact_index, const int time_stage_begin, 
                          const int time_stage_end) {
  contact_sequence_.activateContact(contact_index, time_stage_begin, 
                                    time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void OCP::deactivateContact(const int contact_index, const int time_stage_begin, 
                            const int time_stage_end) {
  contact_sequence_.deactivateContact(contact_index, time_stage_begin, 
                                      time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void OCP::activateContacts(const std::vector<int>& contact_indices, 
                           const int time_stage_begin, 
                           const int time_stage_end) {
  contact_sequence_.activateContacts(contact_indices, time_stage_begin, 
                                     time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void OCP::deactivateContacts(const std::vector<int>& contact_indices, 
                             const int time_stage_begin, 
                             const int time_stage_end) {
  contact_sequence_.deactivateContacts(contact_indices, time_stage_begin, 
                                       time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void OCP::setContactPoint(
    const std::vector<Eigen::Vector3d>& contact_points) {
  assert(contact_points.size() == robots_[0].max_point_contacts());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<robots_.size(); ++i) {
    robots_[i].setContactPoints(contact_points);
  }
}


void OCP::setContactPointByKinematics(const Eigen::VectorXd& q) {
  assert(q.size() == robots_[0].dimq());
  const int dimv = robots_[0].dimv();
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<robots_.size(); ++i) {
    robots_[i].updateKinematics(q, Eigen::VectorXd::Zero(dimv), 
                                Eigen::VectorXd::Zero(dimv));
    robots_[i].setContactPointsByCurrentKinematics();
  }
}


void OCP::clearLineSearchFilter() {
  filter_.clear();
}


double OCP::KKTError() {
  Eigen::VectorXd error = Eigen::VectorXd::Zero(N_+1);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i < N_) {
      error(i) = split_ocps_[i].squaredNormKKTResidual(dtau_);
    }
    else {
      error(N_) = terminal_ocp_.squaredNormKKTResidual();
    }
  }
  return std::sqrt(error.sum());
}


void OCP::computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i == 0) {
      split_ocps_[i].computeKKTResidual(robots_[robot_id], 
                                        contact_sequence_.contactStatus(i), 
                                        t+(i+1)*dtau_, dtau_, q, s_[i], s_[i+1]);
    }
    else if (i < N_) {
      split_ocps_[i].computeKKTResidual(robots_[robot_id], 
                                        contact_sequence_.contactStatus(i), 
                                        t+(i+1)*dtau_, dtau_, s_[i-1].q, s_[i], 
                                        s_[i+1]);
    }
    else {
      terminal_ocp_.computeKKTResidual(robots_[robot_id], t+T_, s_[N_]);
    }
  }
}


void OCP::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
    std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
    std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
    std::cout << "f[" << i << "] = ";
    for (int j=0; j<s_[i].f.size(); ++j) {
      std::cout << s_[i].f[j].transpose() << "; ";
    }
    std::cout << std::endl;
    std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    std::cout << "mu[" << i << "] = " << s_[i].mu_stack().transpose() << std::endl;
  }
  std::cout << "q[" << N_ << "] = " << s_[N_].q.transpose() << std::endl;
  std::cout << "v[" << N_ << "] = " << s_[N_].v.transpose() << std::endl;
}


bool OCP::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    const bool feasible = split_ocps_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time step " << i << std::endl;
      return false;
    }
  }
  return true;
}


void OCP::initConstraints() {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].initConstraints(robots_[robot_id], i, dtau_, s_[i]);
  }
}

} // namespace idocp