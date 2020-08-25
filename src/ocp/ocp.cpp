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
    violations_(Eigen::VectorXd::Zero(N)),
    contact_sequence_(N, std::vector<bool>(robot.max_point_contacts(), false)) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  for (int i=0; i<=N; ++i) {
    robot.normalizeConfiguration(s_[i].q);
  }
  bool feasible = isCurrentSolutionFeasible();
  initConstraints();
}


OCP::OCP() 
  : split_ocps_(),
    terminal_ocp_(),
    robots_(),
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
    violations_(),
    contact_sequence_() {
}


OCP::~OCP() {
}


void OCP::updateSolution(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v, const bool use_line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    if (i < N_) {
      const int robot_id = omp_get_thread_num();
      robots_[robot_id].setContactStatus(contact_sequence_[i]);
      split_ocps_[i].linearizeOCP(robots_[robot_id], t+i*dtau_, dtau_, 
                                  s_[i], s_[i+1]);
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
    split_ocps_[i].computeCondensedDirection(robots_[robot_id], dtau_, d_[i]);
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
              = split_ocps_[i].costAndViolation(robots_[robot_id], t+i*dtau_, 
                                                dtau_, s_[i]);
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
          robots_[robot_id].setContactStatus(contact_sequence_[i]);
          const std::pair<double, double> filter_pair
              = split_ocps_[i].costAndViolation(robots_[robot_id], 
                                                primal_step_size, t+i*dtau_, 
                                                dtau_, s_[i], d_[i], 
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


void OCP::getInitialControlInput(Eigen::VectorXd& u) {
  assert(u.size() == robots_[0].dimv());
  u = s_[0].u;
}


void OCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv) {
  // assert(Kq.rows() == robots_[0].dimv());
  // assert(Kq.cols() == robots_[0].dimv());
  // assert(Kv.rows() == robots_[0].dimv());
  // assert(Kv.cols() == robots_[0].dimv());
  // split_ocps_[0].getStateFeedbackGain(Kq, Kv);
}


bool OCP::setStateTrajectory(const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  Eigen::VectorXd q_normalized = q;
  robots_[0].normalizeConfiguration(q_normalized);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].v = v;
    s_[i].q = q_normalized;
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
    s_[i].v = v0 + i * a;
    robots_[0].integrateConfiguration(q0, v, (double)i, s_[i].q);
  }
  bool feasible = isCurrentSolutionFeasible();
  if (feasible) {
    initConstraints();
  }
  return feasible;
}


void OCP::setContactSequence(
    const std::vector<std::vector<bool>>& contact_sequence) {
  if (contact_sequence.size() != N_) {
    std::cout << "invalid number of the contact sequence: it must be " 
              << N_ << std::endl;
    return;
  }
  for (int i=0; i<N_; ++i) {
    if (contact_sequence[i].size() != robots_[0].max_point_contacts()) {
      std::cout << "invalid dimension of the contact sequence at time step" 
                << i << ": it must be " << robots_[0].max_point_contacts() 
                << std::endl;
      return;
    }
  }
  contact_sequence_ = contact_sequence;
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    const int robot_id = omp_get_thread_num();
    robots_[robot_id].setContactStatus(contact_sequence_[i]);
    s_[i].setContactStatus(robots_[robot_id]);
  }
}


void OCP::setContactPoint(
    const std::vector<Eigen::Vector3d>& contact_points) {
  assert(contact_points.size() == robots_[0].max_point_contacts());
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<contact_points.size(); ++i) {
    robots_[i].setContactPoints(contact_points);
  }
}


void OCP::clearLineSearchFilter() {
  filter_.clear();
}


double OCP::KKTError(const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  Eigen::VectorXd error = Eigen::VectorXd::Zero(N_+1);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N_; ++i) {
    const int robot_id = omp_get_thread_num();
    if (i < N_) {
      robots_[robot_id].setContactStatus(contact_sequence_[i]);
      error(i) = split_ocps_[i].squaredKKTErrorNorm(robots_[robot_id], 
                                                    t+(i+1)*dtau_, dtau_, 
                                                    s_[i], s_[i+1]);
    }
    else {
      error(N_) = terminal_ocp_.squaredKKTErrorNorm(robots_[robot_id], t+T_, 
                                                    s_[N_]);
    }
  }
  return std::sqrt(error.sum());
}


void OCP::printSolution() const {
  for (int i=0; i<N_; ++i) {
    std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
    std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
    std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
    std::cout << "f[" << i << "] = " << s_[i].f.transpose() << std::endl;
    std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    std::cout << "mu[" << i << "] = " << s_[i].mu.transpose() << std::endl;
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
    robots_[robot_id].setContactStatus(contact_sequence_[i]);
    split_ocps_[i].initConstraints(robots_[robot_id], i, dtau_, s_[i]);
  }
}

} // namespace idocp