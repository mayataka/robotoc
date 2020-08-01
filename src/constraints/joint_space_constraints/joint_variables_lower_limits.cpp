#include "constraints/joint_space_constraints/joint_variables_lower_limits.hpp"
#include "constraints/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {

JointVariablesLowerLimits::JointVariablesLowerLimits(const Robot& robot, 
                                                     const Eigen::VectorXd& xmin, 
                                                     const double barrier)
  : dimv_(robot.dimv()),
    dimc_(xmin.size()),
    dim_passive_(robot.dim_passive()),
    has_floating_base_(robot.has_floating_base()),
    barrier_(barrier),
    xmin_(xmin),
    slack_(-xmin_-Eigen::VectorXd::Constant(dimc_, barrier)),
    dual_(Eigen::VectorXd::Constant(dimc_, barrier)),
    residual_(Eigen::VectorXd::Zero(dimc_)),
    duality_(Eigen::VectorXd::Zero(dimc_)),
    dslack_(Eigen::VectorXd::Zero(dimc_)), 
    ddual_(Eigen::VectorXd::Zero(dimc_)) {
  assert(barrier_ > 0);
  assert(xmin_.maxCoeff() < 0);
}


bool JointVariablesLowerLimits::isFeasible(const Eigen::VectorXd& x) {
  assert(x.size() >= dimv_);
  for (int i=0; i<dimc_; ++i) {
    if (x.tail(dimc_).coeff(i) < xmin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVariablesLowerLimits::setSlackAndDual(const double dtau, 
                                                const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  slack_ = dtau * (x.tail(dimc_)-xmin_);
  pdipmfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointVariablesLowerLimits::condenseSlackAndDual(const double dtau, 
                                                     const Eigen::VectorXd& x,
                                                     Eigen::MatrixXd& Cxx, 
                                                     Eigen::VectorXd& Cx) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  assert(Cxx.rows() == dimv_);
  assert(Cxx.cols() == dimv_);
  assert(Cx.size() == dimv_);
  for (int i=0; i<dimv_; ++i) {
    Cxx.coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * dual_.coeff(i) / slack_.coeff(i);
  }
  residual_ = dtau * (xmin_-x.tail(dimc_)) + slack_;
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, duality_);
  Cx.tail(dimc_).array() 
      -= dtau * (dual_.array()*residual_.array()-duality_.array()) 
              / slack_.array();
}


void JointVariablesLowerLimits::computeSlackAndDualDirection(
    const double dtau, const Eigen::VectorXd& dx) {
  assert(dtau > 0);
  assert(dx.size() == dimv_);
  dslack_ = dtau * dx.tail(dimc_) - residual_;
  pdipmfunc::ComputeDualDirection(dual_, slack_, dslack_, duality_, ddual_);
}


double JointVariablesLowerLimits::maxSlackStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


double JointVariablesLowerLimits::maxDualStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


void JointVariablesLowerLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


void JointVariablesLowerLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


double JointVariablesLowerLimits::costSlackBarrier() {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


double JointVariablesLowerLimits::costSlackBarrier(const double step_size) {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


void JointVariablesLowerLimits::augmentDualResidual(const double dtau, 
                                                    Eigen::VectorXd& Cx) {
  assert(dtau > 0);
  assert(Cx.size() == dimv_);
  Cx.tail(dimc_).noalias() -= dtau * dual_;
}


double JointVariablesLowerLimits::residualL1Nrom(const double dtau, 
                                                 const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (xmin_-x.tail(dimc_)) + slack_;
  return residual_.lpNorm<1>();
}


double JointVariablesLowerLimits::residualSquaredNrom(const double dtau, 
                                                      const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (xmin_-x.tail(dimc_)) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, residual_);
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipm
} // namespace idocp