#include "constraints/fb/joint_velocity_lower_limits_fb.hpp"
#include "constraints/fb/fb_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace fb {

JointVelocityLowerLimits::JointVelocityLowerLimits(const Robot& robot, 
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointVelocityLimit().size()),
    barrier_(barrier),
    vmin_(-robot.jointVelocityLimit()),
    slack_(-vmin_-Eigen::VectorXd::Constant(vmin_.size(), barrier)),
    dual_(Eigen::VectorXd::Constant(vmin_.size(), barrier)),
    residual_(Eigen::VectorXd::Zero(vmin_.size())),
    dslack_(Eigen::VectorXd::Zero(vmin_.size())), 
    ddual_(Eigen::VectorXd::Zero(vmin_.size())),
    fb_residual_(Eigen::VectorXd::Zero(vmin_.size())),
    slack_tilde_(Eigen::VectorXd::Zero(vmin_.size())),
    dual_tilde_(Eigen::VectorXd::Zero(vmin_.size())) {
  assert(barrier_ > 0);
  assert(vmin_.maxCoeff() < 0);
}


bool JointVelocityLowerLimits::isFeasible(const Robot& robot, 
                                          const Eigen::VectorXd& v) {
  assert(v.size() == robot.dimv());
  for (int i=0; i<dimc_; ++i) {
    if (v.coeff(i) < vmin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityLowerLimits::setSlackAndDual(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  slack_ = dtau * (v-vmin_);
  fbfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointVelocityLowerLimits::condenseSlackAndDual(const Robot& robot, 
                                                    const double dtau, 
                                                    const Eigen::VectorXd& v,
                                                    Eigen::MatrixXd& Cvv, 
                                                    Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  fbfunc::ComputeFischerBurmeisterRadius(dimc_, barrier_, slack_, dual_, fb_residual_);
  slack_tilde_.array() = 1 - slack_.array() / fb_residual_.array();
  dual_tilde_.array() = 1 - dual_.array() / fb_residual_.array();
  fb_residual_.noalias() -= slack_;
  fb_residual_.noalias() -= dual_;
  for (int i=0; i<dimv_; ++i) {
    Cvv.coeffRef(i, i) += dtau * dtau * slack_tilde_.coeff(i) / dual_tilde_.coeff(i);
  }
  residual_ = dtau * (vmin_-v) + slack_;
  Cv.array() -= dtau * slack_tilde_.array() * residual_.array() / dual_tilde_.array();
  Cv.array() -= dtau * fb_residual_.array() / dual_tilde_.array();
}


void JointVelocityLowerLimits::computeSlackAndDualDirection(
    const Robot& robot, const double dtau, const Eigen::VectorXd& dv) {
  assert(dtau > 0);
  assert(dv.size() == robot.dimv());
  dslack_ = dtau * dv - residual_;
  fbfunc::ComputeDualDirection(dual_tilde_, slack_tilde_, dslack_, fb_residual_, ddual_);
}


double JointVelocityLowerLimits::maxSlackStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return fbfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


double JointVelocityLowerLimits::maxDualStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return fbfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


void JointVelocityLowerLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


void JointVelocityLowerLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


double JointVelocityLowerLimits::costSlackBarrier() {
  return fbfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


double JointVelocityLowerLimits::costSlackBarrier(const double step_size) {
  return fbfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


void JointVelocityLowerLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau, 
                                                   Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(Cv.size() == robot.dimv());
  Cv.noalias() -= dtau * dual_;
}


double JointVelocityLowerLimits::residualL1Nrom(const Robot& robot, 
                                                const double dtau, 
                                                const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (vmin_-v) + slack_;
  return residual_.lpNorm<1>();
}


double JointVelocityLowerLimits::residualSquaredNrom(const Robot& robot, 
                                                     const double dtau, 
                                                     const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (vmin_-v) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  fbfunc::ComputeFischerBurmeisterRadius(dimc_, barrier_, slack_, dual_, residual_);
  residual_.noalias() -= slack_;
  residual_.noalias() -= dual_;
  error += residual_.squaredNorm();
  return error;
}

} // namespace fb
} // namespace idocp