#include "constraints/joint_velocity_lower_limits.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {

JointVelocityLowerLimits::JointVelocityLowerLimits(const Robot& robot)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointVelocityLimit().size()),
    vmin_(-robot.jointVelocityLimit()),
    slack_(-vmin_+Eigen::VectorXd::Constant(vmin_.size(), 1.0e-04)),
    dual_(Eigen::VectorXd::Constant(vmin_.size(), 1.0e-04)),
    residual_(Eigen::VectorXd::Zero(vmin_.size())),
    FB_residual_(Eigen::VectorXd::Zero(vmin_.size())),
    dslack_(Eigen::VectorXd::Zero(vmin_.size())), 
    ddual_(Eigen::VectorXd::Zero(vmin_.size())) {
  for (int i=0; i<vmin_.size(); ++i) {
    assert(vmin_.coeff(i) <= 0);
  }
}


void JointVelocityLowerLimits::setSlackAndDual(const Robot& robot, 
                                               const double barrier,
                                               const Eigen::VectorXd& v) {
  assert(barrier > 0);
  assert(v.size() == robot.dimv());
  slack_ = v - vmin_;
  for (int i=0; i<dimc_; ++i) {
    while (slack_.coeff(i) < barrier) {
      slack_.coeffRef(i) += barrier;
    }
  }
  dual_.array() = (barrier*barrier) / slack_.array();
}


void JointVelocityLowerLimits::linearizeConstraint(const Robot& robot, 
                                                   const double barrier, 
                                                   const Eigen::VectorXd& v) {
  assert(barrier > 0);
  assert(v.size() == robot.dimv());
  residual_ = vmin_ - v;
  for (int i=0; i<dimc_; ++i) {
    const double r = std::sqrt(slack_.coeff(i)*slack_.coeff(i) 
                               +dual_.coeff(i)*dual_.coeff(i) 
                               +2*barrier*barrier);
    FB_residual_.coeffRef(i)  = r - slack_.coeff(i) - dual_.coeff(i);
    dslack_.coeffRef(i) = 1 - slack_.coeff(i) / r;
    ddual_.coeffRef(i) = 1 - dual_.coeff(i) / r;
  }
}


void JointVelocityLowerLimits::updateSlackAndDual(const Robot& robot, 
                                                  const Eigen::VectorXd& dv) {
  assert(dv.size() == robot.dimv());
  residual_.noalias() += slack_ - dv;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointVelocityLowerLimits::condenseSlackAndDual(const Robot& robot, 
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


void JointVelocityLowerLimits::augmentDualResidual(const Robot& robot, 
                                                   Eigen::VectorXd& Cv) {
  assert(Cv.size() == robot.dimv());
  Cv.noalias() -= dual_;
}


double JointVelocityLowerLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const Eigen::VectorXd& v) {
  assert(dtau > 0);
  assert(v.size() == robot.dimv());
  residual_ = vmin_ - v + slack_;
  FB_residual_.array() = slack_.array() * dual_.array();
  double error = 0;
  error += residual_.squaredNorm();
  error += FB_residual_.squaredNorm();
  return error;
}

} // namespace idocp