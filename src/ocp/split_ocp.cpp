#include "idocp/ocp/split_ocp.hpp"

#include <assert.h>


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, 
                   const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    robot_dynamics_(robot),
    riccati_gain_(robot),
    riccati_factorizer_(robot),
    riccati_inverter_(robot),
    Ginv_full_(Eigen::MatrixXd::Zero(
        robot.dimv()+2*robot.max_dimf()+robot.dim_passive(), 
        robot.dimv()+2*robot.max_dimf()+robot.dim_passive())),
    s_tmp_(robot),
    dimv_(robot.dimv()),
    dim_passive_(robot.dim_passive()),
    dimf_(0),
    dimc_(robot.dim_passive()),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.max_point_contacts() > 0) {
    use_kinematics_ = true;
  }
}


SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    robot_dynamics_(),
    riccati_gain_(),
    riccati_factorizer_(),
    riccati_inverter_(),
    Ginv_full_(),
    s_tmp_(),
    dimv_(0),
    dim_passive_(0),
    dimf_(0),
    dimc_(0),
    use_kinematics_(false) {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


void SplitOCP::initConstraints(Robot& robot, const int time_step, 
                               const double dtau, const SplitSolution& s) { 
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void SplitOCP::linearizeOCP(Robot& robot, const ContactStatus& contact_status,  
                            const double t, const double dtau, 
                            const Eigen::VectorXd& q_prev, 
                            const SplitSolution& s, 
                            const SplitSolution& s_next) {
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  setContactStatusForRiccatiRecursion(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // condensing the inverse dynamics
  kkt_residual_.lu.setZero();
  kkt_matrix_.Quu.setZero();
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s.u,
                                    kkt_residual_.lu);
  cost_->luu(robot, cost_data_, t, dtau, s.u, kkt_matrix_.Quu);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s.u, 
                                     kkt_matrix_.Quu, kkt_residual_.lu);
  robot_dynamics_.condenseRobotDynamics(robot, contact_status, dtau, s, 
                                        kkt_matrix_, kkt_residual_);
  // construct the KKT matrix and KKT residual
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  riccati_factorizer_.setStateEquationDerivative(kkt_matrix_.Fqq);
  kkt_matrix_.Qvq() = kkt_matrix_.Qqv().transpose();
  if (contact_status.hasActiveContacts()) {
    kkt_matrix_.Qfa() = kkt_matrix_.Qaf().transpose();
  }
}


void SplitOCP::backwardRiccatiRecursion(
    const double dtau, const RiccatiFactorization& riccati_next, 
    RiccatiFactorization& riccati) {
  assert(dtau > 0);
  riccati_factorizer_.factorize_F(dtau, riccati_next.Pqq, riccati_next.Pqv, 
                                  riccati_next.Pvq, riccati_next.Pvv, 
                                  kkt_matrix_.Qqq(), kkt_matrix_.Qqv(), 
                                  kkt_matrix_.Qvq(), kkt_matrix_.Qvv());
  riccati_factorizer_.factorize_H(dtau, riccati_next.Pqv, riccati_next.Pvv, 
                                  kkt_matrix_.Qaq().transpose(), 
                                  kkt_matrix_.Qav().transpose());
  riccati_factorizer_.factorize_G(dtau, riccati_next.Pvv, kkt_matrix_.Qaa());
  riccati_factorizer_.factorize_la(dtau, riccati_next.Pvq, riccati_next.Pvv, 
                                   kkt_residual_.Fq(), kkt_residual_.Fv(), 
                                   riccati_next.sv, kkt_residual_.la());
  // Computes the matrix inversion
  riccati_inverter_.invert(kkt_matrix_.Qafaf(), kkt_matrix_.Caf(), Ginv_());
  // Computes the state feedback gain and feedforward terms
  riccati_gain_.computeFeedbackGain(Ginv_(), kkt_matrix_.Qafqv(), 
                                    kkt_matrix_.Cqv());
  riccati_gain_.computeFeedforward(Ginv_(), kkt_residual_.laf(), 
                                   kkt_residual_.C());
  // Computes the Riccati factorization matrices
  // Qaq() means Qqa().transpose(). This holds for Qav(), Qfq(), Qfv().
  riccati.Pqq = kkt_matrix_.Qqq();
  riccati.Pqq.noalias() += riccati_gain_.Kaq().transpose() * kkt_matrix_.Qaq();
  riccati.Pqv = kkt_matrix_.Qqv();
  riccati.Pqv.noalias() += riccati_gain_.Kaq().transpose() * kkt_matrix_.Qav();
  riccati.Pvv = kkt_matrix_.Qvv();
  riccati.Pvv.noalias() += riccati_gain_.Kav().transpose() * kkt_matrix_.Qav();
  // Computes the Riccati factorization vectors
  riccati.sq = riccati_next.sq - kkt_residual_.lq();
  riccati.sq.noalias() -= riccati_next.Pqq * kkt_residual_.Fq();
  riccati.sq.noalias() -= riccati_next.Pqv * kkt_residual_.Fv();
  riccati.sq.noalias() -= kkt_matrix_.Qaq().transpose() * riccati_gain_.ka();
  riccati.sv = dtau * riccati_next.sq + riccati_next.sv - kkt_residual_.lv();
  riccati.sv.noalias() -= dtau * riccati_next.Pqq * kkt_residual_.Fq();
  riccati.sv.noalias() -= riccati_next.Pvq * kkt_residual_.Fq();
  riccati.sv.noalias() -= dtau * riccati_next.Pqv * kkt_residual_.Fv();
  riccati.sv.noalias() -= riccati_next.Pvv * kkt_residual_.Fv();
  riccati.sv.noalias() -= kkt_matrix_.Qav().transpose() * riccati_gain_.ka();
  if (dimf_ > 0) {
    riccati.Pqq.noalias() += riccati_gain_.Kfq().transpose() * kkt_matrix_.Qfq();
    riccati.Pqv.noalias() += riccati_gain_.Kfq().transpose() * kkt_matrix_.Qfv();
    riccati.Pvv.noalias() += riccati_gain_.Kfv().transpose() * kkt_matrix_.Qfv();
    riccati.sq.noalias() -= kkt_matrix_.Qfq().transpose() * riccati_gain_.kf();
    riccati.sv.noalias() -= kkt_matrix_.Qfv().transpose() * riccati_gain_.kf();
  }
  if (dimc_ > 0) {
    riccati.Pqq.noalias() += riccati_gain_.Kmuq().transpose() * kkt_matrix_.Cq();
    riccati.Pqv.noalias() += riccati_gain_.Kmuq().transpose() * kkt_matrix_.Cv();
    riccati.Pvv.noalias() += riccati_gain_.Kmuv().transpose() * kkt_matrix_.Cv();
    riccati.sq.noalias() -= kkt_matrix_.Cq().transpose() * riccati_gain_.kmu();
    riccati.sv.noalias() -= kkt_matrix_.Cv().transpose() * riccati_gain_.kmu();
  }
  riccati.Pvq = riccati.Pqv.transpose();
}


void SplitOCP::forwardRiccatiRecursion(const double dtau, SplitDirection& d,   
                                       SplitDirection& d_next) {
  assert(dtau > 0);
  d.da() = riccati_gain_.ka();
  d.da().noalias() += riccati_gain_.Kaq() * d.dq();
  d.da().noalias() += riccati_gain_.Kav() * d.dv();
  d_next.dq() = d.dq() + dtau * d.dv() + kkt_residual_.Fq();
  d_next.dv() = d.dv() + dtau * d.da() + kkt_residual_.Fv();
}


void SplitOCP::computeCondensedDirection(Robot& robot, const double dtau, 
                                         const SplitSolution& s, 
                                         SplitDirection& d) {
  assert(dtau > 0);
  if (dimf_ > 0) {
    d.df() = riccati_gain_.kf();
    d.df().noalias() += riccati_gain_.Kfq() * d.dq();
    d.df().noalias() += riccati_gain_.Kfv() * d.dv();
  }
  if (dimc_ > 0) {
    d.dmu() = riccati_gain_.kmu();
    d.dmu().noalias() += riccati_gain_.Kmuq() * d.dq();
    d.dmu().noalias() += riccati_gain_.Kmuv() * d.dv();
  }
  robot_dynamics_.computeCondensedDirection(dtau, kkt_matrix_, kkt_residual_, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, s, d);
}

 
double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> SplitOCP::costAndConstraintViolation(
    Robot& robot, const double t, const double dtau, const SplitSolution& s) {
  assert(dtau > 0);
  return std::make_pair(cost(robot, t, dtau, s), constraintViolation(dtau));
}


std::pair<double, double> SplitOCP::costAndConstraintViolation(
    Robot& robot, const ContactStatus& contact_status, const double step_size, 
    const double t, const double dtau, const SplitSolution& s, 
    const SplitDirection& d, const SplitSolution& s_next, 
    const SplitDirection& d_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  s_tmp_.setContactStatus(contact_status);
  s_tmp_.a = s.a + step_size * d.da();
  if (contact_status.hasActiveContacts()) {
    s_tmp_.f_stack() = s.f_stack() + step_size * d.df();
    s_tmp_.set_f();
    robot.setContactForces(contact_status, s_tmp_.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp_.q);
  s_tmp_.v = s.v + step_size * d.dv();
  s_tmp_.u = s.u + step_size * d.du;
  if (use_kinematics_) {
    robot.updateKinematics(s_tmp_.q, s_tmp_.v, s_tmp_.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, 
                                             s_tmp_);
  stateequation::ComputeForwardEulerResidual(robot, step_size, dtau, s_tmp_,  
                                             s_next.q, s_next.v, d_next.dq(), 
                                             d_next.dv(), kkt_residual_);
  robot_dynamics_.computeRobotDynamicsResidual(robot, contact_status, dtau, 
                                               s_tmp_, kkt_residual_);
  return std::make_pair(cost(robot, t, dtau, s_tmp_), constraintViolation(dtau));
}


void SplitOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  constraints_->updateDual(constraints_data_, step_size);
}


void SplitOCP::updatePrimal(Robot& robot, const double step_size, 
                            const double dtau,  
                            const RiccatiFactorization& riccati, 
                            const SplitDirection& d, SplitSolution& s) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  s.lmd.noalias() 
      += step_size * (riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq);
  s.gmm.noalias() 
      += step_size * (riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv);
  robot.integrateConfiguration(d.dq(), step_size, s.q);
  s.v.noalias() += step_size * d.dv();
  s.a.noalias() += step_size * d.da();
  s.f_stack().noalias() += step_size * d.df();
  s.set_f();
  s.u.noalias() += step_size * d.du;
  s.beta.noalias() += step_size * d.dbeta;
  s.mu_stack().noalias() += step_size * d.dmu();
  s.set_mu_contact();
  constraints_->updateSlack(constraints_data_, step_size);
}


void SplitOCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                    Eigen::MatrixXd& Kv) const {
  assert(Kq.cols() == dimv_);
  assert(Kq.rows() == dimv_);
  assert(Kv.cols() == dimv_);
  assert(Kv.rows() == dimv_);
  robot_dynamics_.getStateFeedbackGain(riccati_gain_.Kaq(), riccati_gain_.Kav(), 
                                       riccati_gain_.Kfq(), riccati_gain_.Kfv(), 
                                       Kq, Kv);
}


void SplitOCP::computeKKTResidual(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const double t, const double dtau, 
                                  const Eigen::VectorXd& q_prev, 
                                  const SplitSolution& s,
                                  const SplitSolution& s_next) {
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  kkt_residual_.setZeroMinimum();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s.u,
                                    kkt_residual_.lu);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  robot_dynamics_.linearizeRobotDynamics(robot, contact_status, dtau, s, 
                                         kkt_matrix_, kkt_residual_);
}


double SplitOCP::squaredNormKKTResidual(const double dtau) const {
  double error = 0;
  error += kkt_residual_.lq().squaredNorm();
  error += kkt_residual_.lv().squaredNorm();
  error += kkt_residual_.la().squaredNorm();
  error += kkt_residual_.lf().squaredNorm();
  error += kkt_residual_.lu.squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual_);
  error += robot_dynamics_.squaredNormRobotDynamicsResidual(dtau, kkt_residual_);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


double SplitOCP::cost(Robot& robot, const double t, const double dtau, 
                      const SplitSolution& s) {
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  return cost;
}


double SplitOCP::constraintViolation(const double dtau) const {
  double violation = 0;
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual_);
  violation += robot_dynamics_.l1NormRobotDynamicsResidual(dtau, kkt_residual_);
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp