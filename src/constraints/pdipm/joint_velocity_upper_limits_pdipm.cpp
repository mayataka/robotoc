#include "constraints/pdipm/joint_velocity_upper_limits_pdipm.hpp"
#include "constraints/pdipm/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {

JointVelocityUpperLimits::JointVelocityUpperLimits(const Robot& robot, 
                                                   const double barrier) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointVelocityLimit().size()),
    barrier_(barrier),
    vmax_(robot.jointVelocityLimit()),
    slack_(vmax_-Eigen::VectorXd::Constant(vmax_.size(), barrier)),
    dual_(Eigen::VectorXd::Constant(vmax_.size(), barrier)),
    residual_(Eigen::VectorXd::Zero(vmax_.size())),
    duality_(Eigen::VectorXd::Zero(vmax_.size())),
    dslack_(Eigen::VectorXd::Zero(vmax_.size())), 
    ddual_(Eigen::VectorXd::Zero(vmax_.size())) {
  assert(barrier_ > 0);
  assert(vmax_.minCoeff() > 0);
}


bool JointVelocityUpperLimits::isFeasible(const Robot& robot, 
                                          const Eigen::VectorXd& v) {
  assert(v.size() == robot.dimv());
  for (int i=0; i<dimc_; ++i) {
    if (v.coeff(i) > vmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVelocityUpperLimits::setSlackAndDual(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  slack_ = dtau * (vmax_-v);
  pdipmfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointVelocityUpperLimits::condenseSlackAndDual(const Robot& robot, 
                                                    const double dtau, 
                                                    const Eigen::VectorXd& v, 
                                                    Eigen::MatrixXd& Cvv, 
                                                    Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cvv.coeffRef(i, i) += dtau * dtau * dual_.coeff(i) / slack_.coeff(i);
  }
  residual_ = dtau * (v-vmax_) + slack_;
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, duality_);
  Cv.array() += dtau * (dual_.array()*residual_.array()-duality_.array()) / slack_.array();
}


void JointVelocityUpperLimits::computeSlackAndDualDirection(
    const Robot& robot, const double dtau, const Eigen::VectorXd& dv) {
  assert(dtau > 0);
  assert(dv.size() == robot.dimv());
  dslack_ = - dtau * dv - residual_;
  pdipmfunc::ComputeDualDirection(dual_, slack_, dslack_, duality_, ddual_);
}


double JointVelocityUpperLimits::maxSlackStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


double JointVelocityUpperLimits::maxDualStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


void JointVelocityUpperLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


void JointVelocityUpperLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


double JointVelocityUpperLimits::costSlackBarrier() {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


double JointVelocityUpperLimits::costSlackBarrier(const double step_size) {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


void JointVelocityUpperLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau,
                                                   Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(Cv.size() == robot.dimv());
  Cv.noalias() += dtau * dual_;
}


double JointVelocityUpperLimits::residualL1Nrom(const Robot& robot, 
                                                const double dtau, 
                                                const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (v-vmax_) + slack_;
  return residual_.lpNorm<1>();
}


double JointVelocityUpperLimits::residualSquaredNrom(const Robot& robot, 
                                                     const double dtau, 
                                                     const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (v-vmax_) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  residual_.array() = slack_.array() * dual_.array() - barrier_;
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipm
} // namespace idocp