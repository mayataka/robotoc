#include "constraints/soft/joint_velocity_upper_limits_soft.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace softconstraints {

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
    FB_residual_(Eigen::VectorXd::Zero(vmax_.size())),
    dslack_(Eigen::VectorXd::Zero(vmax_.size())), 
    ddual_(Eigen::VectorXd::Zero(vmax_.size())) {
  assert(barrier_ > 0);
  assert(vmax_.array() >= 0);
}


void JointVelocityUpperLimits::setSlackAndDual(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  slack_ = dtau * (vmax_-v);
  for (int i=0; i<dimc_; ++i) {
    while (slack_.coeff(i) < barrier_) {
      slack_.coeffRef(i) += barrier_;
    }
  }
  dual_.array() = (barrier_*barrier_) / slack_.array();
}


void JointVelocityUpperLimits::linearizeConstraint(const Robot& robot, 
                                                   const double dtau, 
                                                   const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (v-vmax_);
  for (int i=0; i<dimc_; ++i) {
    const double r = std::sqrt(slack_.coeff(i)*slack_.coeff(i) 
                               +dual_.coeff(i)*dual_.coeff(i) 
                               +2*barrier_*barrier_);
    FB_residual_.coeffRef(i) = r - slack_.coeff(i) - dual_.coeff(i);
    dslack_.coeffRef(i) = 1 - slack_.coeff(i) / r;
    ddual_.coeffRef(i) = 1 - dual_.coeff(i) / r;
  }
}


void JointVelocityUpperLimits::updateSlackAndDual(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& dv) {
  assert(dtau > 0);
  assert(dv.size() == robot.dimv());
  residual_.noalias() += slack_ + dtau * dv;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointVelocityUpperLimits::condenseSlackAndDual(const Robot& robot, 
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
  Cv.array() += dtau * dslack_.array() * (residual_.array()+slack_.array()) 
                / ddual_.array();
  Cv.array() += dtau * FB_residual_.array() / ddual_.array();
}


void JointVelocityUpperLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau, 
                                                   Eigen::VectorXd& Cv) {
  assert(dtau > 0);
  assert(Cv.size() == robot.dimv());
  // Cv.noalias() += dtau * dual_;
  Cv.noalias() += dual_;
}


double JointVelocityUpperLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = dtau * (v-vmax_) + slack_;
  FB_residual_.array() = slack_.array() * dual_.array();
  double error = 0;
  error += residual_.squaredNorm();
  error += FB_residual_.squaredNorm();
  return error;
}

} // namespace softconstraints
} // namespace idocp