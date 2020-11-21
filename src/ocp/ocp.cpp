#include "idocp/ocp/ocp.hpp"

#include <cmath>
#include <omp.h>
#include <cassert>


namespace idocp {

OCP::OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
         const std::shared_ptr<Constraints>& constraints, const double T, 
         const int N, const int max_num_impulse, const int num_proc)
  : split_ocps_(N, SplitOCP(robot, cost, constraints), max_num_impulse, 
                SplitImpulseOCP(robot, cost->getImpulseCostFunction(), 
                                constraints->getImpulseConstraints())),
    terminal_ocp_(robot, cost, constraints),
    robots_(num_proc, robot),
    filter_(),
    T_(T),
    dtau_(T/N),
    N_(N),
    num_proc_(num_proc),
    s_(N+1, SplitSolution(robot), max_num_impulse, ImpulseSplitSolution(robot)),
    d_(N+1, SplitDirection(robot), max_num_impulse, ImpulseSplitDirection(robot)),
    riccati_factorization_(N+1, RiccatiFactorization(robot), max_num_impulse, RiccatiFactorization(robot)),
    constraint_factorization_(max_num_impulse, 
                              StateConstraintRiccatiFactorization(robot, N, max_num_impulse)),
    constraint_factorizer_(robot, max_num_impulse, num_proc),
    riccati_recursion_(robot, T, N, max_num_impulse, num_proc),
    contact_sequence_(robot, T, N) {
  assert(T > 0);
  assert(N > 0);
  assert(num_proc > 0);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<=N; ++i) {
    robot.normalizeConfiguration(s_[i].q);
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
  ocp_linearizer_.linearizeOCP(split_ocps_, terminal_ocp_, robots_, 
                               contact_sequence_, t, q, v, s_,
                               kkt_matrix_, kkt_residual_);
  riccati_recursion_.backwardRiccatiRecursionTerminal(kkt_matrix_, 
                                                      kkt_residual_, 
                                                      riccati_factorization_);
  riccati_recursion_.backwardRiccatiRecursion(contact_sequence_, kkt_matrix_, 
                                              kkt_residual_, 
                                              riccati_factorization_);
  riccati_recursion_.forwardRiccatiRecursionParallel(contact_sequence_, 
                                                     kkt_matrix_, 
                                                     kkt_residual_);
  riccati_recursion_.forwardRiccatiRecursionSerial(contact_sequence_, 
                                                   kkt_matrix_, kkt_residual_, 
                                                   riccati_factorization_);
  riccati_recursion_.backwardStateConstraintFactorization(contact_sequence_,
                                                          kkt_matrix_, 
                                                          constraint_factorization_);

  ocp_direction_calculator_.computeInitialStateDirection(robots_, s_, q, v, d_);
  constraint_factorizer_.computeLagrangeMultiplierDirection(
      contact_sequence_, riccati_factorization_.impulse, 
      constraint_factorization_, d_[0].dx(), d_.impulse);
  ocp_direction_calculator_.computeDirection(split_ocps_, terminal_ocp_, 
                                             robots_, contact_sequence_, 
                                             riccati_recursion_.getFactorizersHandle(),
                                             riccati_factorization_, 
                                             constraint_factorization_, s_, d_);
  const double primal_step_size = ocp_direction_calculator_.maxPrimalStepSize(contact_sequence_);
  const double dual_step_size = ocp_direction_calculator_.maxDualStepSize(contact_sequence_);
  ocp_solution_integrator_.integrate(split_ocps_, terminal_ocp_, robots_, 
                                     contact_sequence_, kkt_matrix_, kkt_residual_, 
                                     primal_step_size, dual_step_size, d_, s_);
} 


const SplitSolution& OCP::getSolution(const int stage) const {
  assert(stage >= 0);
  assert(stage <= N_);
  return s_[stage];
}


void OCP::getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                               Eigen::MatrixXd& Kv) const {
  assert(stage >= 0);
  assert(stage < N_);
  assert(Kq.rows() == robots_[0].dimv());
  assert(Kq.cols() == robots_[0].dimv());
  assert(Kv.rows() == robots_[0].dimv());
  assert(Kv.cols() == robots_[0].dimv());
  // split_ocps_[stage].getStateFeedbackGain(Kq, Kv);
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


void OCP::setContactStatus(const ContactStatus& contact_status) {
  contact_sequence_.setContactStatusUniformly(contact_status);
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    s_[i].setContactStatus(contact_sequence_.contactStatus(i));
    d_[i].setContactStatus(contact_sequence_.contactStatus(i));
  }
}


ContactSequence& OCP::getContactSequenceHandle() {
  return contact_sequence_;
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
  return ocp_linearizer_.KKTError(split_ocps_, terminal_ocp_, 
                                  contact_sequence_, kkt_residual_);
}


void OCP::computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v) {
  ocp_linearizer_.computeKKTResidual(split_ocps_, terminal_ocp_, robots_, 
                                     contact_sequence_, t, q, v, s_,
                                     kkt_matrix_, kkt_residual_);
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