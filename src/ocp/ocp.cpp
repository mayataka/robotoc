#include "idocp/ocp/ocp.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>


namespace idocp {

OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
         const std::shared_ptr<Constraints>& constraints, const double T, 
         const int N, const int max_num_impulse, const int num_proc)
  : robots_(num_proc, robot),
    contact_sequence_(robot, T, N),
    ocp_linearizer_(T, N, max_num_impulse, num_proc),
    riccati_solver_(robot, T, N, max_num_impulse, num_proc),
    split_ocps_(N, max_num_impulse, robot, cost, constraints),
    kkt_matrix_(N, max_num_impulse, robot),
    kkt_residual_(N, max_num_impulse, robot),
    s_(N, max_num_impulse, robot),
    d_(N, max_num_impulse, robot),
    // filter_(),
    N_(N),
    num_proc_(num_proc),
    T_(T),
    dtau_(T/N) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (num_proc <= 0) {
      throw std::out_of_range("invalid value: num_proc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N; ++i) {
    robot.normalizeConfiguration(s_[i].q);
  }
  for (int i=0; i<max_num_impulse; ++i) {
    robot.normalizeConfiguration(s_.impulse[i].q);
    robot.normalizeConfiguration(s_.aux[i].q);
    robot.normalizeConfiguration(s_.lift[i].q);
  }
  initConstraints();
}


OCP::OCP() {
}


OCP::~OCP() {
}


void OCP::updateSolution(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v, const bool use_line_search) {
  assert(q.size() == robots_[0].dimq());
  assert(v.size() == robots_[0].dimv());
  ocp_linearizer_.linearizeOCP(split_ocps_, robots_, contact_sequence_, 
                               t, q, v, s_, kkt_matrix_, kkt_residual_);
  riccati_solver_.computeNewtonDirection(split_ocps_, robots_, 
                                         contact_sequence_, q, v, s_, d_, 
                                         kkt_matrix_, kkt_residual_);
  const double primal_step_size = riccati_solver_.maxPrimalStepSize();
  const double dual_step_size = riccati_solver_.maxDualStepSize();
  if (use_line_search) {
    // TODO: add filter line search method
  }
  ocp_linearizer_.integrateSolution(split_ocps_, robots_, contact_sequence_, 
                                    kkt_matrix_, kkt_residual_, 
                                    primal_step_size, dual_step_size, d_, s_);
} 


const SplitSolution& OCP::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= N_);
  return s_[stage];
}


void OCP::getStateFeedbackGain(const int time_stage, Eigen::MatrixXd& Kq, 
                               Eigen::MatrixXd& Kv) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  riccati_solver_.getStateFeedbackGain(time_stage, Kq, Kv);
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
  for (int i=0; i<contact_sequence_.totalNumImpulseStages(); ++i) {
    s_.impulse[i].v = v;
    s_.impulse[i].q = q_normalized;
    s_.aux[i].v = v;
    s_.aux[i].q = q_normalized;
  }
  for (int i=0; i<contact_sequence_.totalNumLiftStages(); ++i) {
    s_.lift[i].v = v;
    s_.lift[i].q = q_normalized;
  }
  initConstraints();
  const bool feasible = isCurrentSolutionFeasible();
  return feasible;
}


void OCP::setContactStatusUniformly(const ContactStatus& contact_status) {
  contact_sequence_.setContactStatusUniformly(contact_status);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


void OCP::setDiscreteEvent(const DiscreteEvent& discrete_event) {
  contact_sequence_.setDiscreteEvent(discrete_event);
  for (int i=0; i<=N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
  for (int i=0; i<contact_sequence_.totalNumImpulseStages(); ++i) {
    s_.impulse[i].setImpulseStatus(contact_sequence_.impulseStatus(i));
    d_.impulse[i].setImpulseStatus(contact_sequence_.impulseStatus(i));
    const int stage = contact_sequence_.timeStageAfterImpulse(i);
    s_.aux[i].setContactStatus(contact_sequence_.contactStatus(stage));
    d_.aux[i].setContactStatus(contact_sequence_.contactStatus(stage));
  }
  for (int i=0; i<contact_sequence_.totalNumLiftStages(); ++i) {
    const int stage = contact_sequence_.timeStageAfterLift(i);
    s_.lift[i].setContactStatus(contact_sequence_.contactStatus(stage));
    d_.lift[i].setContactStatus(contact_sequence_.contactStatus(stage));
  }
}


void OCP::shiftImpulse(const int impulse_index, const double impulse_time) {
  contact_sequence_.shiftImpulse(impulse_index, impulse_time);
  for (int i=0; i<=N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
  for (int i=0; i<contact_sequence_.totalNumImpulseStages(); ++i) {
    s_.impulse[i].setImpulseStatus(contact_sequence_.impulseStatus(i));
    d_.impulse[i].setImpulseStatus(contact_sequence_.impulseStatus(i));
    const int stage = contact_sequence_.timeStageAfterImpulse(i);
    s_.aux[i].setContactStatus(contact_sequence_.contactStatus(stage));
    d_.aux[i].setContactStatus(contact_sequence_.contactStatus(stage));
  }
  for (int i=0; i<contact_sequence_.totalNumLiftStages(); ++i) {
    const int stage = contact_sequence_.timeStageAfterLift(i);
    s_.lift[i].setContactStatus(contact_sequence_.contactStatus(stage));
    d_.lift[i].setContactStatus(contact_sequence_.contactStatus(stage));
  }
}


void OCP::shiftLift(const int lift_index, const double lift_time) {
  contact_sequence_.shiftLift(lift_index, lift_time);
  for (int i=0; i<=N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
  for (int i=0; i<contact_sequence_.totalNumImpulseStages(); ++i) {
    s_.impulse[i].setImpulseStatus(contact_sequence_.impulseStatus(i));
    d_.impulse[i].setImpulseStatus(contact_sequence_.impulseStatus(i));
    const int stage = contact_sequence_.timeStageAfterImpulse(i);
    s_.aux[i].setContactStatus(contact_sequence_.contactStatus(stage));
    d_.aux[i].setContactStatus(contact_sequence_.contactStatus(stage));
  }
  for (int i=0; i<contact_sequence_.totalNumLiftStages(); ++i) {
    const int stage = contact_sequence_.timeStageAfterLift(i);
    s_.lift[i].setContactStatus(contact_sequence_.contactStatus(stage));
    d_.lift[i].setContactStatus(contact_sequence_.contactStatus(stage));
  }
}


void OCP::clearLineSearchFilter() {
  // filter_.clear();
}


double OCP::KKTError() {
  return ocp_linearizer_.KKTError(split_ocps_, contact_sequence_, kkt_residual_);
}


void OCP::computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  ocp_linearizer_.computeKKTResidual(split_ocps_, robots_, contact_sequence_, 
                                     t, q, v, s_, kkt_matrix_, kkt_residual_);
}


void OCP::printSolution(const std::string& name, 
                        const std::vector<int> frames) const {
  if (name == "all") {
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
    }
    std::cout << "q[" << N_ << "] = " << s_[N_].q.transpose() << std::endl;
    std::cout << "v[" << N_ << "] = " << s_[N_].v.transpose() << std::endl;
  }
  if (name == "q") {
    for (int i=0; i<=N_; ++i) {
      std::cout << "q[" << i << "] = " << s_[i].q.transpose() << std::endl;
    }
  }
  if (name == "v") {
    for (int i=0; i<=N_; ++i) {
      std::cout << "v[" << i << "] = " << s_[i].v.transpose() << std::endl;
    }
  }
  if (name == "a") {
    for (int i=0; i<N_; ++i) {
      std::cout << "a[" << i << "] = " << s_[i].a.transpose() << std::endl;
    }
  }
  if (name == "f") {
    for (int i=0; i<N_; ++i) {
      std::cout << "f[" << i << "] = ";
      for (int j=0; j<s_[i].f.size(); ++j) {
        std::cout << s_[i].f[j].transpose() << "; ";
      }
    }
  }
  if (name == "u") {
    for (int i=0; i<N_; ++i) {
      std::cout << "u[" << i << "] = " << s_[i].u.transpose() << std::endl;
    }
  }
  if (name == "end-effector") {
    Robot robot = robots_[0];
    for (int i=0; i<N_; ++i) {
      robot.updateFrameKinematics(s_[i].q);
      for (const auto e : frames) {
      std::cout << "ee[" << i << "][" << e << "] = " 
                << robot.framePosition(e).transpose() << std::endl;
      }
    }
  }
}


bool OCP::isCurrentSolutionFeasible() {
  for (int i=0; i<N_; ++i) {
    const bool feasible = split_ocps_[i].isFeasible(robots_[0], s_[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at time stage " << i << std::endl;
      return false;
    }
  }
  const int num_impulse = contact_sequence_.totalNumImpulseStages();
  for (int i=0; i<num_impulse; ++i) {
    const bool feasible = split_ocps_.impulse[i].isFeasible(robots_[0], 
                                                            s_.impulse[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at impulse " << i << std::endl;
      return false;
    }
  }
  for (int i=0; i<num_impulse; ++i) {
    const bool feasible = split_ocps_.aux[i].isFeasible(robots_[0], s_.aux[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at aux " << i << std::endl;
      return false;
    }
  }
  const int num_lift = contact_sequence_.totalNumLiftStages();
  for (int i=0; i<num_lift; ++i) {
    const bool feasible = split_ocps_.lift[i].isFeasible(robots_[0], s_.lift[i]);
    if (!feasible) {
      std::cout << "INFEASIBLE at lift " << i << std::endl;
      return false;
    }
  }
  return true;
}


void OCP::initConstraints() {
  ocp_linearizer_.initConstraints(split_ocps_, robots_, contact_sequence_, s_);
}

} // namespace idocp