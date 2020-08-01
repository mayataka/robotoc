#include "ocp/split_terminal_ocp.hpp"

#include <assert.h>


namespace idocp {

SplitTerminalOCP::SplitTerminalOCP(const Robot& robot, 
                                   const CostFunctionInterface* cost, 
                                   const ConstraintsInterface* constraints) 
  : cost_(const_cast<CostFunctionInterface*>(cost)),
    constraints_(const_cast<ConstraintsInterface*>(constraints)),
    joint_constraints_(robot),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    lq_(Eigen::VectorXd::Zero(robot.dimv())),
    lv_(Eigen::VectorXd::Zero(robot.dimv())),
    q_res_(Eigen::VectorXd::Zero(robot.dimv())),
    v_res_(Eigen::VectorXd::Zero(robot.dimv())),
    // The following variables are only needed for line search
    q_tmp_(Eigen::VectorXd::Zero(robot.dimq())), 
    v_tmp_(Eigen::VectorXd::Zero(robot.dimv())) {
}


SplitTerminalOCP::~SplitTerminalOCP() {
}


bool SplitTerminalOCP::isFeasible(Robot& robot, const Eigen::VectorXd& q, 
                                  const Eigen::VectorXd& v) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  // return joint_constraints_.isFeasible(robot, q, v, a, u);
  return true;
}


void SplitTerminalOCP::initConstraints(Robot& robot, const int time_step,
                                       const double dtau, 
                                       const Eigen::VectorXd& q, 
                                       const Eigen::VectorXd& v) {
  assert(time_step >= 0);
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(a.size() == robot.dimv());
  assert(u.size() == robot.dimv());
  // joint_constraints_.setTimeStep(time_step);
  // joint_constraints_.setSlackAndDual(robot, dtau, q, v, a, u);
}


void SplitTerminalOCP::linearizeOCP(Robot& robot, const double t, 
                                    const Eigen::VectorXd& lmd, 
                                    const Eigen::VectorXd& gmm, 
                                    const Eigen::VectorXd& q, 
                                    const Eigen::VectorXd& v, 
                                    Eigen::MatrixXd& Qqq, Eigen::MatrixXd& Qqv, 
                                    Eigen::MatrixXd& Qvq, Eigen::MatrixXd& Qvv, 
                                    Eigen::VectorXd& Qq, Eigen::VectorXd& Qv) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(Qqq.rows() == robot.dimv());
  assert(Qqq.cols() == robot.dimv());
  assert(Qqv.rows() == robot.dimv());
  assert(Qqv.cols() == robot.dimv());
  assert(Qvq.rows() == robot.dimv());
  assert(Qvq.cols() == robot.dimv());
  assert(Qvv.rows() == robot.dimv());
  assert(Qvv.cols() == robot.dimv());
  assert(Qq.size() == robot.dimv());
  assert(Qv.size() == robot.dimv());
  if (robot.has_floating_base()) {
    cost_->setConfigurationJacobian(robot, q);
  }
  cost_->phiq(t, q, v, lq_);
  cost_->phiv(t, q, v, lv_);
  Qq = - lq_ + lmd;
  Qv = - lv_ + gmm;
  cost_->phiqq(t, q, v, Qqq);
  cost_->phivv(t, q, v, Qvv);
}


void SplitTerminalOCP::computeCondensedDirection(Robot& robot, 
                                                 const double dtau, 
                                                 const Eigen::VectorXd& dq, 
                                                 const Eigen::VectorXd& dv) {
  joint_constraints_.computeSlackAndDualDirection(robot, dtau, dq, dv, da_, du_);
}

 
double SplitTerminalOCP::maxPrimalStepSize() {
  return joint_constraints_.maxSlackStepSize();
}


double SplitTerminalOCP::maxDualStepSize() {
  return joint_constraints_.maxDualStepSize();
}


double SplitTerminalOCP::costDerivativeDotDirection(Robot& robot, 
                                                    const double t, 
                                                    const Eigen::VectorXd& q, 
                                                    const Eigen::VectorXd& v, 
                                                    const Eigen::VectorXd& dq,
                                                    const Eigen::VectorXd& dv) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  if (robot.has_floating_base()) {
    cost_->setConfigurationJacobian(robot, q);
  }
  cost_->phiq(t, q, v, lq_);
  cost_->phiv(t, q, v, lv_);
  double product = 0;
  product += lq_.dot(dq);
  product += lv_.dot(dv);
  return product;
}


double SplitTerminalOCP::terminalCost(Robot& robot, const double t, 
                                      const Eigen::VectorXd& q, 
                                      const Eigen::VectorXd& v) {
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  return cost_->phi(t, q, v);
}


double SplitTerminalOCP::terminalCost(Robot& robot, const double step_size, 
                                      const double t, const Eigen::VectorXd& q, 
                                      const Eigen::VectorXd& v, 
                                      const Eigen::VectorXd& dq, 
                                      const Eigen::VectorXd& dv) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  q_tmp_ = q + step_size * dq;
  v_tmp_ = v + step_size * dv;
  return cost_->phi(t, q_tmp_, v_tmp_);
}


void SplitTerminalOCP::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  joint_constraints_.updateDual(step_size);
}


void SplitTerminalOCP::updatePrimal(Robot& robot, const double step_size, 
                                    const Eigen::VectorXd& dq, 
                                    const Eigen::VectorXd& dv, 
                                    const Eigen::MatrixXd& Pqq, 
                                    const Eigen::MatrixXd& Pqv, 
                                    const Eigen::MatrixXd& Pvq, 
                                    const Eigen::MatrixXd& Pvv, 
                                    const Eigen::VectorXd& sq, 
                                    const Eigen::VectorXd& sv, 
                                    Eigen::VectorXd& q, Eigen::VectorXd& v, 
                                    Eigen::VectorXd& lmd, 
                                    Eigen::VectorXd& gmm) const {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(Pqq.rows() == robot.dimv());
  assert(Pqq.cols() == robot.dimv());
  assert(Pqv.rows() == robot.dimv());
  assert(Pqv.cols() == robot.dimv());
  assert(Pvq.rows() == robot.dimv());
  assert(Pvq.cols() == robot.dimv());
  assert(Pvv.rows() == robot.dimv());
  assert(Pvv.cols() == robot.dimv());
  assert(sq.size() == robot.dimv());
  assert(sv.size() == robot.dimv());
  assert(q.rows() == robot.dimq());
  assert(v.rows() == robot.dimv());
  assert(lmd.rows() == robot.dimv());
  assert(gmm.rows() == robot.dimv());
  q.noalias() += step_size * dq;
  v.noalias() += step_size * dv;
  lmd.noalias() += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm.noalias() += step_size * (Pvq * dq + Pvv * dv - sv);
}


double SplitTerminalOCP::squaredKKTErrorNorm(Robot& robot, const double t, 
                                             const Eigen::VectorXd& lmd, 
                                             const Eigen::VectorXd& gmm, 
                                             const Eigen::VectorXd& q, 
                                             const Eigen::VectorXd& v) {
  assert(lmd.size() == robot.dimv());
  assert(gmm.size() == robot.dimv());
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  // Compute the partial derivatives of the Lagrangian with respect to the 
  // terminal configuration and velocity.
  if (robot.has_floating_base()) {
    cost_->setConfigurationJacobian(robot, q);
  }
  cost_->phiq(t, q, v, lq_);
  cost_->phiv(t, q, v, lv_);
  lq_.noalias() -= lmd;
  lv_.noalias() -= gmm;
  double error = 0;
  error += lq_.squaredNorm();
  error += lv_.squaredNorm();
  return error;
}

} // namespace idocp