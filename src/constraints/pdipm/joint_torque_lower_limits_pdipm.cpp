#include "constraints/pdipm/joint_torque_lower_limits_pdipm.hpp"
#include "constraints/pdipm/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {

JointTorqueLowerLimits::JointTorqueLowerLimits(const Robot& robot, 
                                               const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointEffortLimit().size()),
    barrier_(barrier),
    umin_(-robot.jointEffortLimit()),
    slack_(-umin_-Eigen::VectorXd::Constant(umin_.size(), barrier)),
    dual_(Eigen::VectorXd::Constant(umin_.size(), barrier)),
    residual_(Eigen::VectorXd::Zero(umin_.size())),
    dslack_(Eigen::VectorXd::Zero(umin_.size())), 
    ddual_(Eigen::VectorXd::Zero(umin_.size())) {
  assert(barrier_ > 0);
  for (int i=0; i<umin_.size(); ++i) {
    assert(umin_(i) <= 0);
  }
}


void JointTorqueLowerLimits::setSlackAndDual(const Robot& robot, 
                                             const double dtau, 
                                             const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  slack_ = dtau * (u-umin_);
  pdipmfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointTorqueLowerLimits::condenseSlackAndDual(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& u,
                                                  Eigen::MatrixXd& Cuu, 
                                                  Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cuu.rows() == robot.dimv());
  assert(Cuu.cols() == robot.dimv());
  assert(Cu.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cuu.coeffRef(i, i) += dtau * dtau * dual_.coeff(i) / slack_.coeff(i);
  }
  residual_ = dtau * (umin_-u);
  Cu.array() += dtau * (dual_.array()*residual_.array()+barrier_) 
                / slack_.array();
}


double JointTorqueLowerLimits::computeDirectionAndMaxStepSize(
    const Robot& robot, const double dtau, const Eigen::VectorXd& du) {
  dslack_ = - slack_ + dtau * du - residual_;
  pdipmfunc::ComputeDualDirection(barrier_, dual_, slack_, dslack_, ddual_);
  const double step_size_slack = pdipmfunc::FractionToBoundary(dimc_, slack_, 
                                                               dslack_);
  const double step_size_dual = pdipmfunc::FractionToBoundary(dimc_, dual_, 
                                                              ddual_);
  return std::min(step_size_slack, step_size_dual);
}


void JointTorqueLowerLimits::updateSlackAndDual(const Robot& robot, 
                                                const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
  dual_.noalias() += step_size * ddual_;
}


void JointTorqueLowerLimits::augmentDualResidual(const Robot& robot, 
                                                 const double dtau, 
                                                 Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cu.size() == robot.dimv());
  Cu.noalias() += dtau * dual_;
}


double JointTorqueLowerLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  residual_ = dtau * (umin_-u) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  residual_.array() = slack_.array() * dual_.array();
  error += residual_.squaredNorm();
  return error;
}



} // namespace pdipm
} // namespace idocp