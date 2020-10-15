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
    contact_sequence_(robot, N),
    filter_(),
    T_(T),
    dtau_(T/N),
    step_size_reduction_rate_(0.75),
    min_step_size_(0.05),
    N_(N),
    num_proc_(num_proc),
    s_(N, SplitSolution(robot)),
    s_new_(N, SplitSolution(robot)),
    d_(N, SplitDirection(robot)),
    aux_mat_old_(N, Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    primal_step_sizes_(Eigen::VectorXd::Zero(N)),
    dual_step_sizes_(Eigen::VectorXd::Zero(N)),
    costs_(Eigen::VectorXd::Zero(N)), 
    violations_(Eigen::VectorXd::Zero(N)) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N; ++i) {
    const int robot_id = omp_get_thread_num();
    robot.normalizeConfiguration(s_[i].q);
    robot.normalizeConfiguration(s_new_[i].q);
  }
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
    d_(),
    aux_mat_old_(),
    primal_step_sizes_(),
    dual_step_sizes_(),
    costs_(), 
    violations_(),
    contact_sequence_() {
}


ParNMPC::~ParNMPC() {
}


void ParNMPC::updateSolution(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, 
                             const bool use_line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i == 0) {
      split_ocps_[i].coarseUpdate(robots_[robot_id], 
                                  contact_sequence_.contactStatus(i), 
                                  t+(i+1)*dtau_, dtau_, q, v, s_[i], s_[i+1], 
                                  aux_mat_old_[i+1], d_[i], s_new_[i]);
    }
    else if (i < N_-1) {
      split_ocps_[i].coarseUpdate(robots_[robot_id], 
                                  contact_sequence_.contactStatus(i), 
                                  t+(i+1)*dtau_, dtau_, s_[i-1].q, s_[i-1].v, 
                                  s_[i], s_[i+1], aux_mat_old_[i+1], d_[i], 
                                  s_new_[i]);
    }
    else {
      split_ocps_[i].coarseUpdateTerminal(robots_[robot_id], 
                                          contact_sequence_.contactStatus(i), 
                                          t+(i+1)*dtau_, dtau_, s_[i-1].q, 
                                          s_[i-1].v, s_[i], d_[i], s_new_[i]);
    }
  }
  for (int i=N_-2; i>=0; --i) {
    split_ocps_[i].backwardCorrectionSerial(robots_[0], s_[i+1], s_new_[i+1], 
                                            s_new_[i]);
  }
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=N_-2; i>=0; --i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].backwardCorrectionParallel(robots_[robot_id], d_[i], 
                                              s_new_[i]);
  }
  for (int i=1; i<N_; ++i) {
    split_ocps_[i].forwardCorrectionSerial(robots_[0], s_[i-1], s_new_[i-1], 
                                           s_new_[i]);
  }
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i > 0) {
      split_ocps_[i].forwardCorrectionParallel(robots_[robot_id], d_[i], 
                                               s_new_[i]);
    }
    split_ocps_[i].computePrimalAndDualDirection(robots_[robot_id], dtau_, 
                                                 s_[i], s_new_[i], d_[i]);
    primal_step_sizes_.coeffRef(i) = split_ocps_[i].maxPrimalStepSize();
    dual_step_sizes_.coeffRef(i) = split_ocps_[i].maxDualStepSize();
  }
  double primal_step_size = primal_step_sizes_.minCoeff();
  const double dual_step_size = dual_step_sizes_.minCoeff();
  if (use_line_search) {
    if (filter_.isEmpty()) {
      #pragma omp parallel for num_threads(num_proc_)
      for (int i=0; i<N_; ++i) {
        const int robot_id = omp_get_thread_num();
        if (i < N_-1) {
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndConstraintViolation(
                    robots_[robot_id], t+(i+1)*dtau_, dtau_, s_[i]);
          costs_.coeffRef(i) = filter_pair.first;
          violations_.coeffRef(i) = filter_pair.second;
        }
        else {
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndConstraintViolationTerminal(
                  robots_[robot_id], t+(i+1)*dtau_, dtau_, s_[i]);
          costs_.coeffRef(i) = filter_pair.first;
          violations_.coeffRef(i) = filter_pair.second;
        }
      } 
      filter_.augment(costs_.sum(), violations_.sum());
    }
    while (primal_step_size > min_step_size_) {
      #pragma omp parallel for num_threads(num_proc_)
      for (int i=0; i<N_; ++i) {
        const int robot_id = omp_get_thread_num();
        if (i == 0) {
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndConstraintViolation(
                  robots_[robot_id], contact_sequence_.contactStatus(i), 
                  primal_step_size, t+(i+1)*dtau_, dtau_, q, v, s_[i], d_[i], 
                  s_new_[i]);
          costs_.coeffRef(i) = filter_pair.first;
          violations_.coeffRef(i) = filter_pair.second;
        }
        else if (i < N_-1) {
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndConstraintViolation(
                  robots_[robot_id], contact_sequence_.contactStatus(i), 
                  primal_step_size, t+(i+1)*dtau_, dtau_, s_[i-1], d_[i-1], 
                  s_[i], d_[i], s_new_[i]);
          costs_.coeffRef(i) = filter_pair.first;
          violations_.coeffRef(i) = filter_pair.second;
        }
        else {
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndConstraintViolationTerminal(
                  robots_[robot_id], contact_sequence_.contactStatus(i), 
                  primal_step_size, t+(i+1)*dtau_, dtau_, s_[i-1], d_[i-1], 
                  s_[i], d_[i], s_new_[i]);
          costs_.coeffRef(i) = filter_pair.first;
          violations_.coeffRef(i) = filter_pair.second;
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
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].updatePrimal(robots_[robot_id], primal_step_size, dtau_,
                                d_[i], s_[i]);
    split_ocps_[i].updateDual(dual_step_size);
    split_ocps_[i].getAuxiliaryMatrix(aux_mat_old_[i]);
  }
} 



void ParNMPC::getControlInput(const int stage, Eigen::VectorXd& u) const {
  assert(stage >= 0);
  assert(stage < N_);
  assert(u.size() == robots_[0].dimv());
  u = s_[stage].u;
}


void ParNMPC::getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                                   Eigen::MatrixXd& Kv) const {
  assert(stage >= 0);
  assert(stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  split_ocps_[stage].getStateFeedbackGain(Kq, Kv);
}


bool ParNMPC::setStateTrajectory(const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  Eigen::VectorXd q_normalized = q;
  robots_[0].normalizeConfiguration(q_normalized);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].q = q_normalized;
    s_[i].v = v;
    s_new_[i].q = q_normalized;
    s_new_[i].v = v;
  }
  initConstraints();
  const bool feasible = isCurrentSolutionFeasible();
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
    s_[i].v = v0 + i * a;
    s_new_[i].v = s_[i].v;
    robots_[0].integrateConfiguration(q0, v, (double)i, s_[i].q);
    s_new_[i].q = s_[i].q;
  }
  initConstraints();
  const bool feasible = isCurrentSolutionFeasible();
  return feasible;
}


void ParNMPC::setAuxiliaryMatrixGuessByTerminalCost(const double t) {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].getTerminalCostHessian(robots_[robot_id], t+(i+1)*dtau_, 
                                          s_[i], aux_mat_old_[i]);
  }
}


void ParNMPC::activateContact(const int contact_index, 
                              const int time_stage_begin, 
                              const int time_stage_end) {
  contact_sequence_.activateContact(contact_index, time_stage_begin, 
                                    time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    s_new_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void ParNMPC::deactivateContact(const int contact_index, 
                                const int time_stage_begin, 
                                const int time_stage_end) {
  contact_sequence_.deactivateContact(contact_index, time_stage_begin, 
                                      time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    s_new_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void ParNMPC::activateContacts(const std::vector<int>& contact_indices, 
                               const int time_stage_begin, 
                               const int time_stage_end) {
  contact_sequence_.activateContacts(contact_indices, time_stage_begin, 
                                     time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    s_new_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void ParNMPC::deactivateContacts(const std::vector<int>& contact_indices, 
                                 const int time_stage_begin, 
                                 const int time_stage_end) {
  contact_sequence_.deactivateContacts(contact_indices, time_stage_begin, 
                                       time_stage_end);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    s_new_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void ParNMPC::setContactPoint(
    const std::vector<Eigen::Vector3d>& contact_points) {
  assert(contact_points.size() == robots_[0].max_point_contacts());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<robots_.size(); ++i) {
    robots_[i].setContactPoints(contact_points);
  }
}


void ParNMPC::setContactPointByKinematics(const Eigen::VectorXd& q) {
  assert(q.size() == robots_[0].dimq());
  const int dimv = robots_[0].dimv();
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<robots_.size(); ++i) {
    robots_[i].updateKinematics(q, Eigen::VectorXd::Zero(dimv), 
                                Eigen::VectorXd::Zero(dimv));
    robots_[i].setContactPointsByCurrentKinematics();
  }
}


void ParNMPC::clearLineSearchFilter() {
  filter_.clear();
}


double ParNMPC::KKTError() {
  Eigen::VectorXd error = Eigen::VectorXd::Zero(N_);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    error(i) = split_ocps_[i].squaredNormKKTResidual(dtau_);
  }
  return std::sqrt(error.sum());
}


void ParNMPC::computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                                const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i == 0) {
      split_ocps_[i].computeKKTResidual(robots_[robot_id], 
                                        contact_sequence_.contactStatus(i), 
                                        t+(i+1)*dtau_, dtau_, q, v, s_[i], 
                                        s_[i+1]);
    }
    else if (i<N_-1) {
      split_ocps_[i].computeKKTResidual(robots_[robot_id], 
                                        contact_sequence_.contactStatus(i), 
                                        t+(i+1)*dtau_, dtau_, s_[i-1].q, 
                                        s_[i-1].v, s_[i], s_[i+1]);
    }
    else {
      split_ocps_[i].computeKKTResidualTerminal(
          robots_[robot_id], contact_sequence_.contactStatus(i), t+(i+1)*dtau_, 
          dtau_, s_[i-1].q, s_[i-1].v, s_[i]);
    }
  }
}


void ParNMPC::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
    std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
    std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
    std::cout << "f[" << i << "] = ";
    for (int j=0; j<s_[i].f.size(); ++j) {
      std::cout << s_[i].f[j].transpose() << "; ";
    }
    std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    std::cout << "mu[" << i << "] = " << s_[i].mu_stack().transpose() << std::endl;
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
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    split_ocps_[i].initConstraints(robots_[robot_id], i+1, dtau_, s_[i]);
  }
}

} // namespace idocp