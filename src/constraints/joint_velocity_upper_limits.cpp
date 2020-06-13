#include "constraints/joint_velocity_upper_limits.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {

JointVelocityUpperLimits::JointVelocityUpperLimits(const Robot& robot)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointVelocityLimit().size()),
    vmax_(robot.jointVelocityLimit()),
    slack_(vmax_-Eigen::VectorXd::Constant(vmax_.size(), 1.0e-04)),
    dual_(Eigen::VectorXd::Constant(vmax_.size(), 1.0e-04)),
    residual_(Eigen::VectorXd::Zero(vmax_.size())),
    FB_residual_(Eigen::VectorXd::Zero(vmax_.size())),
    dslack_(Eigen::VectorXd::Zero(vmax_.size())), 
    ddual_(Eigen::VectorXd::Zero(vmax_.size())) {
  for (int i=0; i<vmax_.size(); ++i) {
    assert(vmax_.coeff(i) >= 0);
  }
}


void JointVelocityUpperLimits::setSlackAndDual(const Robot& robot, 
                                               const double barrier,
                                               const Eigen::VectorXd& v) {
  assert(barrier > 0);
  assert(v.size() == robot.dimv());
  slack_ = vmax_ - v;
  for (int i=0; i<dimc_; ++i) {
    while (slack_.coeff(i) < barrier) {
      slack_.coeffRef(i) += barrier;
    }
  }
  dual_.array() = (barrier*barrier) / slack_.array();
}


void JointVelocityUpperLimits::linearizeConstraint(const Robot& robot, 
                                                   const double barrier, 
                                                   const Eigen::VectorXd& v) {
  assert(barrier > 0);
  assert(v.size() == robot.dimv());
  residual_ = v - vmax_;
  for (int i=0; i<dimc_; ++i) {
    const double r = std::sqrt(slack_.coeff(i)*slack_.coeff(i) 
                               +dual_.coeff(i)*dual_.coeff(i) 
                               +2*barrier*barrier);
    FB_residual_.coeffRef(i)  = r - slack_.coeff(i) - dual_.coeff(i);
    dslack_.coeffRef(i) = 1 - slack_.coeff(i) / r;
    ddual_.coeffRef(i) = 1 - dual_.coeff(i) / r;
  }
}


void JointVelocityUpperLimits::updateSlackAndDual(const Robot& robot, 
                                                  const Eigen::VectorXd& dv) {
  assert(dv.size() == robot.dimv());
  residual_.noalias() += slack_ + dv;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointVelocityUpperLimits::condenseSlackAndDual(const Robot& robot, 
                                                    Eigen::MatrixXd& Cvv, 
                                                    Eigen::VectorXd& Cv) const {
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cvv.coeffRef(i, i) += dslack_.coeff(i) / ddual_.coeff(i);
  }
  Cv.array() += dslack_.array() * (residual_.array()+slack_.array()) 
                / ddual_.array();
  Cv.array() += FB_residual_.array() / ddual_.array();
}


void JointVelocityUpperLimits::augmentDualResidual(const Robot& robot, 
                                                   Eigen::VectorXd& Cv) {
  assert(Cv.size() == robot.dimv());
  Cv.noalias() += dual_;
}


double JointVelocityUpperLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const Eigen::VectorXd& v) {
  assert(v.size() == robot.dimv());
  residual_ = v - vmax_ + slack_;
  FB_residual_.array() = slack_.array() * dual_.array();
  double error = 0;
  error += residual_.squaredNorm();
  error += FB_residual_.squaredNorm();
  return error;
}

} // namespace idocp