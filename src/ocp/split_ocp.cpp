#include "idocp/ocp/split_ocp.hpp"

#include <assert.h>


namespace idocp {

SplitOCP::SplitOCP(const Robot& robot, 
                   const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(robot),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    kkt_residual_(robot),
    kkt_matrix_(robot),
    state_equation_(robot),
    robot_dynamics_(robot),
    riccati_gain_(robot),
    riccati_factorizer_(robot),
    riccati_inverter_(robot),
    Ginv_(Eigen::MatrixXd::Zero(
        robot.dimv()+2*robot.max_dimf()+robot.dim_passive(), 
        robot.dimv()+2*robot.max_dimf()+robot.dim_passive())),
    s_tmp_(robot),
    dimv_(robot.dimv()),
    dimf_(robot.dimf()),
    dimc_(robot.dim_passive()+robot.dimf()) {
}


SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    kkt_residual_(),
    kkt_matrix_(),
    state_equation_(),
    robot_dynamics_(),
    riccati_gain_(),
    riccati_factorizer_(),
    riccati_inverter_(),
    Ginv_(),
    s_tmp_(),
    dimv_(0),
    dimf_(0),
    dimc_(0) {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(const Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


void SplitOCP::initConstraints(const Robot& robot, const int time_step, 
                               const double dtau, const SplitSolution& s) { 
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


void SplitOCP::linearizeOCP(Robot& robot, const double t, const double dtau, 
                            const SplitSolution& s, 
                            const SplitSolution& s_next) {
  assert(dtau > 0);
  setContactStatus(robot);
  kkt_residual_.setContactStatus(robot);
  kkt_matrix_.setContactStatus(robot);
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  // condensing the inverse dynamics
  kkt_residual_.lu.setZero();
  kkt_matrix_.Quu.setZero();
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_.lu);
  cost_->luu(robot, cost_data_, t, dtau, s.u, kkt_matrix_.Quu);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s.u, 
                                     kkt_matrix_.Quu, kkt_residual_.lu);
  robot_dynamics_.condenseRobotDynamics(robot, dtau, s, kkt_matrix_, 
                                        kkt_residual_);
  // forms the KKT matrix and KKT residual
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  state_equation_.linearizeForwardEuler(robot, dtau, s, s_next, kkt_matrix_, 
                                        kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  riccati_factorizer_.setIntegrationSensitivities(kkt_matrix_.Fqq, 
                                                  kkt_matrix_.Fqv);
  riccati_gain_.setContactStatus(robot);
  riccati_inverter_.setContactStatus(robot);
  kkt_matrix_.Qvq() = kkt_matrix_.Qqv().transpose();
  if (robot.has_active_contacts()) {
    kkt_matrix_.Qfa() = kkt_matrix_.Qaf().transpose();
  }
}


void SplitOCP::backwardRiccatiRecursion(
    const double dtau, const RiccatiFactorization& riccati_next, 
    RiccatiFactorization& riccati) {
  assert(dtau > 0);
  riccati_factorizer_.factorizeF(dtau, riccati_next.Pqq, riccati_next.Pqv, 
                                 riccati_next.Pvq, riccati_next.Pvv, 
                                 kkt_matrix_.Qqq(), kkt_matrix_.Qqv(), 
                                 kkt_matrix_.Qvq(), kkt_matrix_.Qvv());
  riccati_factorizer_.factorizeH(dtau, riccati_next.Pqv, riccati_next.Pvv, 
                                 kkt_matrix_.Qaq().transpose(), 
                                 kkt_matrix_.Qav().transpose());
  riccati_factorizer_.factorizeG(dtau, riccati_next.Pvv, kkt_matrix_.Qaa());
  kkt_residual_.la().noalias() += dtau * riccati_next.Pvq * kkt_residual_.Fq(); 
  kkt_residual_.la().noalias() += dtau * riccati_next.Pvv * kkt_residual_.Fv();
  kkt_residual_.la().noalias() -= dtau * riccati_next.sv;
  // Computes the matrix inversion
  riccati_inverter_.invert(kkt_matrix_.Qafaf(), kkt_matrix_.Caf(), 
                           Ginv_active());
  // Computes the state feedback gain and feedforward terms
  riccati_gain_.computeFeedbackGain(Ginv_active(), kkt_matrix_.Qafqv(), 
                                    kkt_matrix_.Cqv());
  riccati_gain_.computeFeedforward(Ginv_active(), kkt_residual_.laf(), 
                                   kkt_residual_.C());
  // Computes the Riccati factorization matrices
  // Qaq() means Qqa().transpose(). This holds for Qav(), Qfq(), Qfv().
  riccati.Pqq = kkt_matrix_.Qqq();
  riccati.Pqq.noalias() += riccati_gain_.Kaq().transpose() * kkt_matrix_.Qaq();
  riccati.Pqv = kkt_matrix_.Qqv();
  riccati.Pqv.noalias() += riccati_gain_.Kaq().transpose() * kkt_matrix_.Qav();
  riccati.Pvq = kkt_matrix_.Qvq();
  riccati.Pvq.noalias() += riccati_gain_.Kav().transpose() * kkt_matrix_.Qaq();
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
    riccati.Pvq.noalias() += riccati_gain_.Kfv().transpose() * kkt_matrix_.Qfq();
    riccati.Pvv.noalias() += riccati_gain_.Kfv().transpose() * kkt_matrix_.Qfv();
    riccati.sq.noalias() -= kkt_matrix_.Qfq().transpose() * riccati_gain_.kf();
    riccati.sv.noalias() -= kkt_matrix_.Qfv().transpose() * riccati_gain_.kf();
  }
  if (dimc_ > 0) {
    riccati.Pqq.noalias() += riccati_gain_.Kmuq().transpose() * kkt_matrix_.Cq();
    riccati.Pqv.noalias() += riccati_gain_.Kmuq().transpose() * kkt_matrix_.Cv();
    riccati.Pvq.noalias() += riccati_gain_.Kmuv().transpose() * kkt_matrix_.Cq();
    riccati.Pvv.noalias() += riccati_gain_.Kmuv().transpose() * kkt_matrix_.Cv();
    riccati.sq.noalias() -= kkt_matrix_.Cq().transpose() * riccati_gain_.kmu();
    riccati.sv.noalias() -= kkt_matrix_.Cv().transpose() * riccati_gain_.kmu();
  }
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


void SplitOCP::computeCondensedDirection(const Robot& robot, const double dtau, 
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
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, d);
}

 
double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


std::pair<double, double> SplitOCP::costAndViolation(Robot& robot, 
                                                     const double t, 
                                                     const double dtau, 
                                                     const SplitSolution& s) {
  assert(dtau > 0);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  cost += constraints_->costSlackBarrier(constraints_data_);
  double violation = 0;
  violation += state_equation_.violationL1Norm(kkt_residual_);
  violation += robot_dynamics_.violationL1Norm(robot, dtau, s, kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, s);
  return std::make_pair(cost, violation);
}


std::pair<double, double> SplitOCP::costAndViolation(
    Robot& robot, const double step_size, const double t, const double dtau, 
    const SplitSolution& s, const SplitDirection& d, 
    const SplitSolution& s_next, const SplitDirection& d_next) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dtau > 0);
  s_tmp_.a = s.a + step_size * d.da();
  if (robot.has_active_contacts()) {
    s_tmp_.f_active() = s.f_active() + step_size * d.df();
    robot.setContactForces(s_tmp_.f);
  }
  robot.integrateConfiguration(s.q, d.dq(), step_size, s_tmp_.q);
  s_tmp_.v = s.v + step_size * d.dv();
  s_tmp_.u = s.u + step_size * d.du;
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s_tmp_.q, s_tmp_.v, s_tmp_.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s_tmp_);
  cost += constraints_->costSlackBarrier(constraints_data_, step_size);
  double violation = 0;
  violation += state_equation_.computeForwardEulerViolationL1Norm(
      robot, step_size, dtau, s_tmp_, s_next.q, s_next.v, 
      d_next.dq(), d_next.dv(), kkt_residual_);
  violation += robot_dynamics_.computeViolationL1Norm(robot, dtau, s_tmp_, 
                                                      kkt_residual_);
  violation += constraints_->residualL1Nrom(robot, constraints_data_, dtau, 
                                            s_tmp_);
  return std::make_pair(cost, violation);
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
  s.f_active().noalias() += step_size * d.df();
  s.u.noalias() += step_size * d.du;
  s.beta.noalias() += step_size * d.dbeta;
  s.mu_active().noalias() += step_size * d.dmu();
  constraints_->updateSlack(constraints_data_, step_size);
}


void SplitOCP::getStateFeedbackGain(Eigen::MatrixXd& Kq, 
                                    Eigen::MatrixXd& Kv) const {
  assert(Kq.cols() == dimv_);
  assert(Kq.rows() == dimv_);
  assert(Kv.cols() == dimv_);
  assert(Kv.rows() == dimv_);
  // Kq = du_dq_ + du_da_ * Kaq_ + du_df_.leftCols(dimf_) * Kfq_.topRows(dimf_);
  // Kv = du_dv_ + du_da_ * Kav_ + du_df_.leftCols(dimf_) * Kfv_.topRows(dimf_);
}


double SplitOCP::squaredKKTErrorNorm(Robot& robot, const double t, 
                                     const double dtau, const SplitSolution& s,
                                     const SplitSolution& s_next) {
  assert(dtau > 0);
  kkt_matrix_.setContactStatus(robot);
  kkt_residual_.setContactStatus(robot);
  kkt_residual_.setZeroMinimum();
  if (robot.has_active_contacts()) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  cost_->lu(robot, cost_data_, t, dtau, s.u, kkt_residual_.lu);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, 
                                    kkt_residual_.lu);
  state_equation_.linearizeForwardEuler(robot, dtau, s, s_next, kkt_matrix_, 
                                        kkt_residual_);
  robot_dynamics_.augmentRobotDynamics(robot, dtau, s, kkt_matrix_, 
                                       kkt_residual_);
  double error = kkt_residual_.squaredKKTErrorNorm(dtau);
  error += constraints_->squaredKKTErrorNorm(robot, constraints_data_, dtau, s);
  return error;
}

} // namespace idocp