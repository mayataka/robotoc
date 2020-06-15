#include "constraints/soft/joint_torque_upper_limits_soft.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace softconstraints {

JointTorqueUpperLimits::JointTorqueUpperLimits(const Robot& robot,
                                               const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointEffortLimit().size()),
    barrier_(barrier),
    umax_(robot.jointEffortLimit()),
    slack_(umax_-Eigen::VectorXd::Constant(umax_.size(), 1.0e-04)),
    dual_(Eigen::VectorXd::Constant(umax_.size(), 1.0e-04)),
    residual_(Eigen::VectorXd::Zero(umax_.size())),
    FB_residual_(Eigen::VectorXd::Zero(umax_.size())),
    dslack_(Eigen::VectorXd::Zero(umax_.size())), 
    ddual_(Eigen::VectorXd::Zero(umax_.size())) {
  assert(barrier_ > 0);
  assert(umax_.array() >= 0);
}


void JointTorqueUpperLimits::setSlackAndDual(const Robot& robot, 
                                             const double dtau, 
                                             const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  slack_ = dtau * (umax_-u);
  for (int i=0; i<dimc_; ++i) {
    while (slack_.coeff(i) < barrier_) {
      slack_.coeffRef(i) += barrier_;
    }
  }
  dual_.array() = (barrier_*barrier_) / slack_.array();
}


void JointTorqueUpperLimits::linearizeConstraint(const Robot& robot, 
                                                 const double dtau, 
                                                 const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  residual_ = dtau * (u-umax_);
  for (int i=0; i<dimc_; ++i) {
    const double r = std::sqrt(slack_.coeff(i)*slack_.coeff(i) 
                               +dual_.coeff(i)*dual_.coeff(i) 
                               +2*barrier_*barrier_);
    FB_residual_.coeffRef(i) = r - slack_.coeff(i) - dual_.coeff(i);
    dslack_.coeffRef(i) = 1 - slack_.coeff(i) / r;
    ddual_.coeffRef(i) = 1 - dual_.coeff(i) / r;
  }
}


void JointTorqueUpperLimits::updateSlackAndDual(const Robot& robot, 
                                                const double dtau, 
                                                const Eigen::VectorXd& du) {
  assert(dtau > 0);
  assert(du.size() == robot.dimv());
  residual_.noalias() += slack_ + dtau * du;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointTorqueUpperLimits::condenseSlackAndDual(const Robot& robot, 
                                                  const double dtau, 
                                                  Eigen::MatrixXd& Cuu, 
                                                  Eigen::VectorXd& Cu) const {
  assert(dtau > 0);
  assert(Cuu.rows() == robot.dimv());
  assert(Cuu.cols() == robot.dimv());
  assert(Cu.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cuu.coeffRef(i, i) += dtau * dtau * dslack_.coeff(i) / ddual_.coeff(i);
  }
  Cu.array() += dtau * dslack_.array() * (residual_.array()+slack_.array()) 
                / ddual_.array();
  Cu.array() += dtau * FB_residual_.array() / ddual_.array();
}


void JointTorqueUpperLimits::augmentDualResidual(const Robot& robot, 
                                                 const double dtau, 
                                                 Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cu.size() == robot.dimv());
  // Cu.noalias() += dtau * dual_;
  Cu.noalias() += dual_;
}


double JointTorqueUpperLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  residual_ = dtau * (u-umax_) + slack_;
  FB_residual_.array() = slack_.array() * dual_.array();
  double error = 0;
  error += residual_.squaredNorm();
  error += FB_residual_.squaredNorm();
  return error;
}

} // namespace softconstraints
} // namespace idocp