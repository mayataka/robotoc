#include "constraints/joint_position_upper_limits.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {

JointPositionUpperLimits::JointPositionUpperLimits(const Robot& robot, 
                                                   const double barrier)
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.upperJointPositionLimit().size()),
    barrier_(barrier),
    qmax_(robot.upperJointPositionLimit()),
    slack_(qmax_-Eigen::VectorXd::Constant(qmax_.size(), barrier)),
    dual_(Eigen::VectorXd::Constant(qmax_.size(), barrier)),
    residual_(Eigen::VectorXd::Zero(qmax_.size())),
    FB_residual_(Eigen::VectorXd::Zero(qmax_.size())),
    dslack_(Eigen::VectorXd::Zero(qmax_.size())), 
    ddual_(Eigen::VectorXd::Zero(qmax_.size())) {
  assert(barrier_ > 0);
}


void JointPositionUpperLimits::setSlackAndDual(const Robot& robot, 
                                               const double dtau, 
                                               const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  slack_ = dtau * (qmax_-q);
  for (int i=0; i<dimc_; ++i) {
    while (slack_.coeff(i) < barrier_) {
      slack_.coeffRef(i) += barrier_;
    }
  }
  dual_.array() = (barrier_*barrier_) / slack_.array();
}


void JointPositionUpperLimits::linearizeConstraint(const Robot& robot, 
                                                   const double dtau, 
                                                   const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  residual_ = dtau * (q-qmax_);
  for (int i=0; i<dimc_; ++i) {
    const double r = std::sqrt(slack_.coeff(i)*slack_.coeff(i) 
                               +dual_.coeff(i)*dual_.coeff(i) 
                               +2*barrier_*barrier_);
    FB_residual_.coeffRef(i) = r - slack_.coeff(i) - dual_.coeff(i);
    dslack_.coeffRef(i) = 1 - slack_.coeff(i) / r;
    ddual_.coeffRef(i) = 1 - dual_.coeff(i) / r;
  }
}


void JointPositionUpperLimits::updateSlackAndDual(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& dq) {
  assert(dtau > 0);
  assert(dq.size() == robot.dimv());
  residual_.noalias() += slack_ + dtau * dq;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointPositionUpperLimits::condenseSlackAndDual(const Robot& robot, 
                                                    const double dtau, 
                                                    Eigen::MatrixXd& Cqq, 
                                                    Eigen::VectorXd& Cq) const {
  assert(dtau > 0);
  assert(Cqq.rows() == robot.dimv());
  assert(Cqq.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cqq.coeffRef(i, i) += dtau * dtau * dslack_.coeff(i) / ddual_.coeff(i);
  }
  Cq.array() += dtau * dslack_.array() * (residual_.array()+slack_.array()) 
                / ddual_.array();
  Cq.array() += dtau * FB_residual_.array() / ddual_.array();
}


void JointPositionUpperLimits::augmentDualResidual(const Robot& robot, 
                                                   const double dtau, 
                                                   Eigen::VectorXd& Cq) {
  assert(dtau > 0);
  assert(Cq.size() == robot.dimv());
  Cq.noalias() += dtau * dual_;
}


double JointPositionUpperLimits::squaredConstraintsResidualNrom(
    const Robot& robot, const double dtau, const Eigen::VectorXd& q) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  residual_ = dtau * (q-qmax_) + slack_;
  FB_residual_.array() = slack_.array() * dual_.array();
  double error = 0;
  error += residual_.squaredNorm();
  error += FB_residual_.squaredNorm();
  return error;
}

} // namespace idocp