#include "constraints/joint_torque_upper_limits.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {

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
    ddual_(Eigen::VectorXd::Zero(umax_.size())),
    newton_residual_(Eigen::VectorXd::Zero(umax_.size())),
    partial_du_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
  assert(barrier_ > 0);
  for (int i=0; i<umax_.size(); ++i) {
    assert(umax_.coeff(i) >= 0);
  }
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
                                                const Eigen::MatrixXd& du_dq,
                                                const Eigen::MatrixXd& du_dv,
                                                const Eigen::MatrixXd& du_da,
                                                const Eigen::VectorXd& dq, 
                                                const Eigen::VectorXd& dv,
                                                const Eigen::VectorXd& da) {
  assert(dtau > 0);
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(dq.size() == robot.dimv());
  assert(dv.size() == robot.dimv());
  assert(da.size() == robot.dimv());
  residual_.noalias() += slack_;
  residual_.noalias() += du_dq * dq;
  residual_.noalias() += du_dv * dv;
  residual_.noalias() += du_da * da;
  ddual_.array() = (dslack_.array()*residual_.array()+FB_residual_.array()) 
                    / (ddual_.array());
  slack_.noalias() -= residual_;
  dual_.noalias() += ddual_;
}


void JointTorqueUpperLimits::condenseSlackAndDual(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::MatrixXd& du_dq,
                                                  const Eigen::MatrixXd& du_dv,
                                                  const Eigen::MatrixXd& du_da, 
                                                  Eigen::MatrixXd& Cqq, 
                                                  Eigen::MatrixXd& Cqv, 
                                                  Eigen::MatrixXd& Cqa, 
                                                  Eigen::MatrixXd& Cvv, 
                                                  Eigen::MatrixXd& Cva, 
                                                  Eigen::MatrixXd& Caa,
                                                  Eigen::VectorXd& Cq, 
                                                  Eigen::VectorXd& Cv, 
                                                  Eigen::VectorXd& Ca) {
  assert(dtau > 0);
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(Cqq.rows() == robot.dimv());
  assert(Cqq.cols() == robot.dimv());
  assert(Cqv.rows() == robot.dimv());
  assert(Cqv.cols() == robot.dimv());
  assert(Cqa.rows() == robot.dimv());
  assert(Cqa.cols() == robot.dimv());
  assert(Cvv.rows() == robot.dimv());
  assert(Cvv.cols() == robot.dimv());
  assert(Cva.rows() == robot.dimv());
  assert(Cva.cols() == robot.dimv());
  assert(Caa.rows() == robot.dimv());
  assert(Caa.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  assert(Ca.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    partial_du_.row(i) = (dtau*dtau*dslack_.coeff(i)/ddual_.coeff(i)) * du_dq.row(i);
  }
  Cqq.noalias() += du_dq.transpose() * partial_du_;
  for (int i=0; i<dimv_; ++i) {
    partial_du_.row(i) = (dtau*dtau*dslack_.coeff(i)/ddual_.coeff(i)) * du_dv.row(i);
  }
  Cqv.noalias() += du_dq.transpose() * partial_du_;
  Cvv.noalias() += du_dv.transpose() * partial_du_;
  for (int i=0; i<dimv_; ++i) {
    partial_du_.row(i) = (dtau*dtau*dslack_.coeff(i)/ddual_.coeff(i)) * du_da.row(i);
  }
  Cqa.noalias() += du_dq.transpose() * partial_du_;
  Cva.noalias() += du_dv.transpose() * partial_du_;
  Caa.noalias() += du_da.transpose() * partial_du_;
  newton_residual_.array() = dtau * dslack_.array() * (residual_.array()+slack_.array()) 
                             / ddual_.array();
  newton_residual_.array() += dtau * FB_residual_.array() / ddual_.array();
  Cq.noalias() += du_dq.transpose() * newton_residual_;
  Cv.noalias() += du_dv.transpose() * newton_residual_;
  Ca.noalias() += du_da.transpose() * newton_residual_;
}


void JointTorqueUpperLimits::augmentDualResidual(const Robot& robot, 
                                                 const double dtau, 
                                                 const Eigen::MatrixXd& du_dq,
                                                 const Eigen::MatrixXd& du_dv, 
                                                 const Eigen::MatrixXd& du_da, 
                                                 Eigen::VectorXd& Cq, 
                                                 Eigen::VectorXd& Cv, 
                                                 Eigen::VectorXd& Ca) {
  assert(dtau > 0);
  assert(du_dq.rows() == robot.dimv());
  assert(du_dq.cols() == robot.dimv());
  assert(du_dv.rows() == robot.dimv());
  assert(du_dv.cols() == robot.dimv());
  assert(du_da.rows() == robot.dimv());
  assert(du_da.cols() == robot.dimv());
  assert(Cq.size() == robot.dimv());
  assert(Cv.size() == robot.dimv());
  assert(Ca.size() == robot.dimv());
  Cq.noalias() += dtau * du_dq.transpose() * dual_;
  Cv.noalias() += dtau * du_dv.transpose() * dual_;
  Ca.noalias() += dtau * du_da.transpose() * dual_;
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

} // namespace idocp