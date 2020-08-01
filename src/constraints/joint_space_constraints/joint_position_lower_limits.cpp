#include "constraints/joint_space_constraints/joint_position_lower_limits.hpp"
#include "constraints/pdipm_func.hpp"

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
    duality_(Eigen::VectorXd::Zero(qmin_.size())),
    dslack_(Eigen::VectorXd::Zero(qmin_.size())), 
    ddual_(Eigen::VectorXd::Zero(qmin_.size())) {
  assert(barrier_ > 0);
}


bool JointPositionLowerLimits::isFeasible(const Robot& robot, 
                                          const Eigen::VectorXd& q) {
  assert(q.size() == robot.dimq());
  for (int i=0; i<dimc_; ++i) {
    if (q.coeff(i) < qmin_.coeff(i)) {
      return false;
    }
  }
  return true;
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
  residual_ = dtau * (qmin_-q) + slack_;
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, duality_);
  Cq.array() -= dtau * (dual_.array()*residual_.array()-duality_.array()) / slack_.array();
}


void JointPositionLowerLimits::computeSlackAndDualDirection(
    const Robot& robot, const double dtau, const Eigen::VectorXd& dq) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  dslack_ = dtau * dq - residual_;
  pdipmfunc::ComputeDualDirection(dual_, slack_, dslack_, duality_, ddual_);
}


double JointPositionLowerLimits::maxSlackStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


double JointPositionLowerLimits::maxDualStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


void JointPositionLowerLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


void JointPositionLowerLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


double JointPositionLowerLimits::costSlackBarrier() {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


double JointPositionLowerLimits::costSlackBarrier(const double step_size) {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


void JointPositionLowerLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau, 
                                                   Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(Cq.size() == robot.dimv());
  Cq.noalias() -= dtau * dual_;
}


double JointPositionLowerLimits::residualL1Nrom(const Robot& robot, 
                                                const double dtau, 
                                                const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  residual_ = dtau * (qmin_-q) + slack_;
  return residual_.lpNorm<1>();
}


double JointPositionLowerLimits::residualSquaredNrom(const Robot& robot, 
                                                     const double dtau, 
                                                     const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  residual_ = dtau * (qmin_-q) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, residual_);
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipm
} // namespace idocp