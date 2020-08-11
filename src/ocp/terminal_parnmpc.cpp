#include "idocp/ocp/terminal_parnmpc.hpp"

#include <assert.h>


namespace idocp {

TerminalParNMPC::TerminalParNMPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost)
  : cost_(cost),
    cost_data_(robot),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    lq_(Eigen::VectorXd::Zero(robot.dimv())),
    lv_(Eigen::VectorXd::Zero(robot.dimv())) {
}


TerminalParNMPC::TerminalParNMPC() 
  : cost_(),
    dimq_(0),
    dimv_(0),
    lq_(),
    lv_() {
}


TerminalParNMPC::~TerminalParNMPC() {
}


void TerminalParNMPC::linearizeOCP(Robot& robot, const double t, 
                                   const Eigen::VectorXd& q, 
                                   const Eigen::VectorXd& v, 
                                   Eigen::VectorXd& phiq, 
                                   Eigen::VectorXd& phiv) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  cost_->phiq(robot, cost_data_, t, q, v, phiq);
  cost_->phiv(robot, cost_data_, t, q, v, phiv);
}


double TerminalParNMPC::terminalCost(Robot& robot, const double t, 
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v) {
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  return cost_->phi(robot, cost_data_, t, q, v);
}


double TerminalParNMPC::terminalCost(Robot& robot, const double step_size, 
                                 const double t, const Eigen::VectorXd& q, 
                                 const Eigen::VectorXd& v, 
                                 const Eigen::VectorXd& dq, 
                                 const Eigen::VectorXd& dv) {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  if (robot.has_floating_base()) {
    q_tmp_ = q;
    robot.integrateConfiguration(dq, step_size, q_tmp_);
  }
  else {
    q_tmp_ = q + step_size * dq;
  }
  v_tmp_ = v + step_size * dv;
  return cost_->phi(robot, cost_data_, t, q_tmp_, v_tmp_);
}


void TerminalParNMPC::updateDual(const double step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  // TODO: add inequality constraints at the terminal OCP.
  // joint_constraints_.updateDual(step_size);
}


void TerminalParNMPC::updatePrimal(Robot& robot, const double step_size, 
                               const Eigen::MatrixXd& Pqq, 
                               const Eigen::MatrixXd& Pqv, 
                               const Eigen::MatrixXd& Pvq, 
                               const Eigen::MatrixXd& Pvv, 
                               const Eigen::VectorXd& sq, 
                               const Eigen::VectorXd& sv, 
                               const Eigen::VectorXd& dq, 
                               const Eigen::VectorXd& dv, Eigen::VectorXd& lmd, 
                               Eigen::VectorXd& gmm, Eigen::VectorXd& q, 
                               Eigen::VectorXd& v) const {
  assert(step_size > 0);
  assert(step_size <= 1);
  assert(Pqq.rows() == dimv_);
  assert(Pqq.cols() == dimv_);
  assert(Pqv.rows() == dimv_);
  assert(Pqv.cols() == dimv_);
  assert(Pvq.rows() == dimv_);
  assert(Pvq.cols() == dimv_);
  assert(Pvv.rows() == dimv_);
  assert(Pvv.cols() == dimv_);
  assert(sq.size() == dimv_);
  assert(sv.size() == dimv_);
  assert(dq.size() == dimv_);
  assert(dv.size() == dimv_);
  assert(lmd.rows() == dimv_);
  assert(gmm.rows() == dimv_);
  assert(q.rows() == dimq_);
  assert(v.rows() == dimv_);
  lmd.noalias() += step_size * (Pqq * dq + Pqv * dv - sq);
  gmm.noalias() += step_size * (Pvq * dq + Pvv * dv - sv);
  robot.integrateConfiguration(dq, step_size, q);
  v.noalias() += step_size * dv;
}


double TerminalParNMPC::squaredKKTErrorNorm(Robot& robot, const double t, 
                                        const Eigen::VectorXd& lmd, 
                                        const Eigen::VectorXd& gmm, 
                                        const Eigen::VectorXd& q, 
                                        const Eigen::VectorXd& v) {
  assert(lmd.size() == dimv_);
  assert(gmm.size() == dimv_);
  assert(q.size() == dimq_);
  assert(v.size() == dimv_);
  // Compute the partial derivatives of the Lagrangian with respect to the 
  // terminal configuration and velocity.
  if (robot.has_floating_base()) {
    robot.computeConfigurationJacobian(q);
  }
  cost_->phiq(robot, cost_data_, t, q, v, lq_);
  cost_->phiv(robot, cost_data_, t, q, v, lv_);
  cost_->phiq(robot, cost_data_, t, q, v, lq_);
  cost_->phiv(robot, cost_data_, t, q, v, lv_);
  lq_.noalias() -= lmd;
  lv_.noalias() -= gmm;
  double error = 0;
  error += lq_.squaredNorm();
  error += lv_.squaredNorm();
  return error;
}

} // namespace idocp