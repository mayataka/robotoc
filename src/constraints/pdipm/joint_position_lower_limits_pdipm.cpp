#include "constraints/pdipm/joint_position_lower_limits_pdipm.hpp"
#include "constraints/pdipm/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {

JointPositionLowerLimits::JointPositionLowerLimits(const Robot& robot, 
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.lowerJointPositionLimit().size()),
    barrier_(barrier),
    qmin_(robot.lowerJointPositionLimit()),
    slack_(-qmin_-Eigen::VectorXd::Constant(qmin_.size(), barrier)),
    dual_(Eigen::VectorXd::Constant(qmin_.size(), barrier)),
    residual_(Eigen::VectorXd::Zero(qmin_.size())),
    dslack_(Eigen::VectorXd::Zero(qmin_.size())), 
    ddual_(Eigen::VectorXd::Zero(qmin_.size())) {
  assert(barrier_ > 0);
}


void JointPositionLowerLimits::setSlackAndDual(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  slack_ = dtau * (q-qmin_);
  pdipmfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointPositionLowerLimits::condenseSlackAndDual(const Robot& robot, 
                                                    const double dtau, 
                                                    const Eigen::VectorXd& q,
                                                    Eigen::MatrixXd& Cqq, 
                                                    Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(Cqq.rows() == robot.dimv());
  assert(Cqq.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cqq.coeffRef(i, i) += dtau * dtau * dual_.coeff(i) / slack_.coeff(i);
  }
  residual_ = dtau * (qmin_-q);
  Cq.array() -= dtau * (dual_.array()*residual_.array()+barrier_) 
                / slack_.array();
}


double JointPositionLowerLimits::computeDirectionAndMaxStepSize(
    const Robot& robot, const double dtau, const Eigen::VectorXd& dq) {
  dslack_ = - slack_ + dtau * dq - residual_;
  pdipmfunc::ComputeDualDirection(barrier_, dual_, slack_, dslack_, ddual_);
  const double step_size_slack = pdipmfunc::FractionToBoundary(dimc_, slack_, 
                                                               dslack_);
  const double step_size_dual = pdipmfunc::FractionToBoundary(dimc_, dual_, 
                                                              ddual_);
  return std::min(step_size_slack, step_size_dual);
}


void JointPositionLowerLimits::updateSlackAndDual(const Robot& robot, 
                                                  const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
  dual_.noalias() += step_size * ddual_;
}


void JointPositionLowerLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau, 
                                                   Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(Cq.size() == robot.dimv());
  Cq.noalias() += dtau * dual_;
}


double JointPositionLowerLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  residual_ = dtau * (qmin_-q) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  residual_.array() = slack_.array() * dual_.array();
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipm
} // namespace idocp