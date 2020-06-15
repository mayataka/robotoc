#include "constraints/soft/joint_velocity_lower_limits_soft.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace softconstraints {

JointVelocityLowerLimits::JointVelocityLowerLimits(const Robot& robot,
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointVelocityLimit().size()),
    barrier_(barrier),
    vmin_(-robot.jointVelocityLimit()),
    slack_(-vmin_+Eigen::VectorXd::Constant(vmin_.size(), 1.0e-04)),
    dual_(Eigen::VectorXd::Constant(vmin_.size(), 1.0e-04)),
    residual_(Eigen::VectorXd::Zero(vmin_.size())),
    FB_residual_(Eigen::VectorXd::Zero(vmin_.size())),
    dslack_(Eigen::VectorXd::Zero(vmin_.size())), 
    ddual_(Eigen::VectorXd::Zero(vmin_.size())) {
  assert(barrier_ > 0);
  assert(vmin_.array() <= 0);
}


void JointVelocityLowerLimits::setSlackAndDual(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  slack_ = dtau * (v-vmin_);
  for (int i=0; i<dimc_; ++i) {
    while (slack_.coeff(i) < barrier_) {
      slack_.coeffRef(i) += barrier_;
    }
  }
  dual_.array() = (barrier_*barrier_) / slack_.array();
}


void JointVelocityLowerLimits::linearizeConstraint(const Robot& robot, 
                                                   const double dtau, 
                                                   const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (vmin_-v);
  for (int i=0; i<dimc_; ++i) {
    const double r = std::sqrt(slack_.coeff(i)*slack_.coeff(i) 
                               +dual_.coeff(i)*dual_.coeff(i) 
                               +2*barrier_*barrier_);
    FB_residual_.coeffRef(i) = r - slack_.coeff(i) - dual_.coeff(i);
    dslack_.coeffRef(i) = 1 - slack_.coeff(i) / r;
    ddual_.coeffRef(i) = 1 - dual_.coeff(i) / r;
  }
}


void JointVelocityLowerLimits::updateSlackAndDual(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& dv) {
  assert(dtau > 0);
  assert(dv.size() == robot.dimv());
  residual_.noalias() += slack_ - dtau * dv;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointVelocityLowerLimits::condenseSlackAndDual(const Robot& robot, 
                                                    const double dtau, 
                                                    Eigen::MatrixXd& Cvv, 
                                                    Eigen::VectorXd& Cv) const {
  assert(dtau > 0);
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cvv.coeffRef(i, i) += dtau * dtau * dslack_.coeff(i) / ddual_.coeff(i);
  }
  Cv.array() -= dtau * dslack_.array() * (residual_.array()+slack_.array()) 
                / ddual_.array();
  Cv.array() -= dtau * FB_residual_.array() / ddual_.array();
}


void JointVelocityLowerLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau, 
                                                   Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(Cv.size() == robot.dimv());
  Cv.noalias() -= dual_;
  // Cv.noalias() -= dtau * dual_;
}


double JointVelocityLowerLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (vmin_-v) + slack_;
  FB_residual_.array() = slack_.array() * dual_.array();
  double error = 0;
  error += residual_.squaredNorm();
  error += FB_residual_.squaredNorm();
  return error;
}

} // namespace softconstraints
} // namespace idocp