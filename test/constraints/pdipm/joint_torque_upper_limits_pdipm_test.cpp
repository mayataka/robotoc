#include "constraints/pdipm/joint_torque_upper_limits_pdipm.hpp"
#include "constraints/pdipm/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {

JointTorqueUpperLimits::JointTorqueUpperLimits(const Robot& robot, 
                                               const double barrier) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimc_(robot.jointEffortLimit().size()),
    barrier_(barrier),
    umax_(robot.jointEffortLimit()),
    slack_(umax_-Eigen::VectorXd::Constant(umax_.size(), barrier)),
    dual_(Eigen::VectorXd::Constant(umax_.size(), barrier)),
    residual_(Eigen::VectorXd::Zero(umax_.size())),
    dslack_(Eigen::VectorXd::Zero(umax_.size())), 
    ddual_(Eigen::VectorXd::Zero(umax_.size())) {
  assert(barrier_ > 0);
  for (int i=0; i<umax_.size(); ++i) {
    assert(umax_(i) >= 0);
  }
}


bool JointTorqueUpperLimits::isFeasible(const Robot& robot, 
                                        const Eigen::VectorXd& u) {
  assert(u.size() == robot.dimv());
  for (int i=0; i<dimc_; ++i) {
    if (u.coeff(i) > umax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointTorqueUpperLimits::setSlackAndDual(const Robot& robot, 
                                             const double dtau, 
                                             const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  slack_ = dtau * (umax_-u);
  pdipmfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointTorqueUpperLimits::condenseSlackAndDual(const Robot& robot, 
                                                  const double dtau, 
                                                  const Eigen::VectorXd& u, 
                                                  Eigen::MatrixXd& Cuu, 
                                                  Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cuu.rows() == robot.dimv());
  assert(Cuu.cols() == robot.dimv());
  assert(Cu.size() == robot.dimv());
  for (int i=0; i<dimv_; ++i) {
    Cuu.coeffRef(i, i) += dtau * dtau * dual_.coeff(i) / slack_.coeff(i);
  }
  residual_ = dtau * (u-umax_);
  Cu.array() += dtau * (dual_.array()*residual_.array()+barrier_) 
                / slack_.array();
}


std::pair<double, double> JointTorqueUpperLimits
    ::computeDirectionAndMaxStepSize(const Robot& robot, const double dtau, 
                                   const Eigen::VectorXd& du) {
  dslack_ = - slack_ - dtau * du - residual_;
  pdipmfunc::ComputeDualDirection(barrier_, dual_, slack_, dslack_, ddual_);
  const double step_size_slack = pdipmfunc::FractionToBoundary(dimc_, slack_, 
                                                               dslack_);
  const double step_size_dual = pdipmfunc::FractionToBoundary(dimc_, dual_, 
                                                              ddual_);
  return std::make_pair(step_size_slack, step_size_dual);
}


void JointTorqueUpperLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


void JointTorqueUpperLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


double JointTorqueUpperLimits::slackBarrier() {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


double JointTorqueUpperLimits::slackBarrier(const double step_size) {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


void JointTorqueUpperLimits::augmentDualResidual(const Robot& robot, 
                                                 const double dtau,
                                                 Eigen::VectorXd& Cu) {
  assert(dtau > 0);
  assert(Cu.size() == robot.dimv());
  Cu.noalias() -= dtau * dual_;
}


double JointTorqueUpperLimits::residualL1Nrom(const Robot& robot, 
                                              const double dtau, 
                                              const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  residual_ = dtau * (u-umax_) + slack_;
  return residual_.lpNorm<1>();
}


double JointTorqueUpperLimits::residualSquaredNrom(const Robot& robot, 
                                                   const double dtau, 
                                                   const Eigen::VectorXd& u) {
  assert(dtau > 0);
  assert(u.size() == robot.dimv());
  residual_ = dtau * (u-umax_) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  residual_.array() = slack_.array() * dual_.array() - barrier_;
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipm 
} // namespace idocp