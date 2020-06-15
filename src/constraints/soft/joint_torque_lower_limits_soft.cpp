#include "constraints/soft/joint_torque_lower_limits_soft.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace softconstraints {

JointTorqueLowerLimits::JointTorqueLowerLimits(const Robot& robot,
                                               const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointEffortLimit().size()),
    barrier_(barrier),
    umin_(-robot.jointEffortLimit()),
    slack_(umin_-Eigen::VectorXd::Constant(umin_.size(), 1.0e-04)),
    dual_(Eigen::VectorXd::Constant(umin_.size(), 1.0e-04)),
    residual_(Eigen::VectorXd::Zero(umin_.size())),
    FB_residual_(Eigen::VectorXd::Zero(umin_.size())),
    dslack_(Eigen::VectorXd::Zero(umin_.size())), 
    ddual_(Eigen::VectorXd::Zero(umin_.size())) {
  assert(barrier_ > 0);
  assert(umin_.array() <= 0);
}


void JointTorqueLowerLimits::setSlackAndDual(const Robot& robot, 
                                             const double dtau, 
                                             const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  slack_ = dtau * (u-umin_);
  for (int i=0; i<dimc_; ++i) {
    while (slack_.coeff(i) < barrier_) {
      slack_.coeffRef(i) += barrier_;
    }
  }
  dual_.array() = (barrier_*barrier_) / slack_.array();
}


void JointTorqueLowerLimits::linearizeConstraint(const Robot& robot, 
                                                 const double dtau, 
                                                 const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  residual_ = dtau * (umin_-u);
  for (int i=0; i<dimc_; ++i) {
    const double r = std::sqrt(slack_.coeff(i)*slack_.coeff(i) 
                               +dual_.coeff(i)*dual_.coeff(i) 
                               +2*barrier_*barrier_);
    FB_residual_.coeffRef(i) = r - slack_.coeff(i) - dual_.coeff(i);
    dslack_.coeffRef(i) = 1 - slack_.coeff(i) / r;
    ddual_.coeffRef(i) = 1 - dual_.coeff(i) / r;
  }
}


void JointTorqueLowerLimits::updateSlackAndDual(const Robot& robot, 
                                                const double dtau, 
                                                const Eigen::VectorXd& du) {
  assert(dtau > 0);
  assert(du.size() == robot.dimv());
  residual_.noalias() += slack_ - dtau * du;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointTorqueLowerLimits::condenseSlackAndDual(const Robot& robot, 
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
  Cu.array() -= dtau * dslack_.array() * (residual_.array()+slack_.array()) 
                / ddual_.array();
  Cu.array() -= dtau * FB_residual_.array() / ddual_.array();
}


void JointTorqueLowerLimits::augmentDualResidual(const Robot& robot, 
                                                 const double dtau, 
                                                 Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cu.size() == robot.dimv());
  // Cu.noalias() -= dtau * dual_;
  Cu.noalias() -= dual_;
}


double JointTorqueLowerLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  residual_ = dtau * (umin_-u) + slack_;
  FB_residual_.array() = slack_.array() * dual_.array();
  double error = 0;
  error += residual_.squaredNorm();
  error += FB_residual_.squaredNorm();
  return error;
}

} // namespace softconstraints
} // namespace idocp