#include "idocp/constraints/joint_space_constraints/joint_variables_upper_limits.hpp"
#include "idocp/constraints/pdipm_func.hpp"

#include <cmath>
#include <assert.h>


namespace idocp {
namespace pdipm {

JointVariablesUpperLimits::JointVariablesUpperLimits(
    const Robot& robot, const Eigen::Ref<const Eigen::VectorXd>& xmax, 
    const double barrier)
  : dimv_(robot.dimv()),
    dimc_(xmax.size()),
    dim_passive_(robot.dim_passive()),
    has_floating_base_(robot.has_floating_base()),
    barrier_(barrier),
    xmax_(xmax),
    slack_(xmax_-Eigen::VectorXd::Constant(dimc_, barrier)),
    dual_(Eigen::VectorXd::Constant(dimc_, barrier)),
    residual_(Eigen::VectorXd::Zero(dimc_)),
    duality_(Eigen::VectorXd::Zero(dimc_)),
    dslack_(Eigen::VectorXd::Zero(dimc_)), 
    ddual_(Eigen::VectorXd::Zero(dimc_)) {
  assert(barrier_ > 0);
  assert(xmax_.minCoeff() > 0);
}


JointVariablesUpperLimits::JointVariablesUpperLimits()
  : dimv_(0),
    dimc_(0),
    dim_passive_(0),
    has_floating_base_(false),
    barrier_(0),
    xmax_(),
    slack_(),
    dual_(),
    residual_(),
    duality_(),
    dslack_(), 
    ddual_() {
}


JointVariablesUpperLimits::~JointVariablesUpperLimits() {
}


bool JointVariablesUpperLimits::isFeasible(
    const Eigen::Ref<const Eigen::VectorXd>& x) {
  assert(x.size() >= dimv_);
  for (int i=0; i<dimc_; ++i) {
    if (x.tail(dimc_).coeff(i) > xmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


void JointVariablesUpperLimits::setSlackAndDual(
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  slack_ = dtau * (xmax_-x.tail(dimc_));
  pdipmfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


void JointVariablesUpperLimits::condenseSlackAndDual(
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& x, 
    Eigen::Ref<Eigen::MatrixXd> Cxx, Eigen::Ref<Eigen::VectorXd> Cx) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  assert(Cxx.rows() == dimv_);
  assert(Cxx.cols() == dimv_);
  assert(Cx.size() == dimv_);
  for (int i=0; i<dimc_; ++i) {
    Cxx.coeffRef(dim_passive_+i, dim_passive_+i) 
        += dtau * dtau * dual_.coeff(i) / slack_.coeff(i);
  }
  residual_ = dtau * (x.tail(dimc_)-xmax_) + slack_;
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, duality_);
  Cx.tail(dimc_).array() 
      += dtau * (dual_.array()*residual_.array()-duality_.array()) 
              / slack_.array();
}


void JointVariablesUpperLimits::computeSlackAndDualDirection(
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& dx) {
  assert(dtau > 0);
  assert(dx.size() == dimv_);
  dslack_ = - dtau * dx.tail(dimc_) - residual_;
  pdipmfunc::ComputeDualDirection(dual_, slack_, dslack_, duality_, ddual_);
}


double JointVariablesUpperLimits::maxSlackStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


double JointVariablesUpperLimits::maxDualStepSize(const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


void JointVariablesUpperLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


void JointVariablesUpperLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


double JointVariablesUpperLimits::costSlackBarrier() {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


double JointVariablesUpperLimits::costSlackBarrier(const double step_size) {
  return pdipmfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


void JointVariablesUpperLimits::augmentDualResidual(
    const double dtau, Eigen::Ref<Eigen::VectorXd> Cx) {
  assert(dtau > 0);
  assert(Cx.size() == dimv_);
  Cx.tail(dimc_).noalias() += dtau * dual_;
}


double JointVariablesUpperLimits::residualL1Nrom(
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (x.tail(dimc_)-xmax_) + slack_;
  return residual_.lpNorm<1>();
}


double JointVariablesUpperLimits::residualSquaredNrom(
    const double dtau, const Eigen::Ref<const Eigen::VectorXd>& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (x.tail(dimc_)-xmax_) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  pdipmfunc::ComputeDualityResidual(barrier_, slack_, dual_, residual_);
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipm
} // namespace idocp