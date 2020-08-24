#ifndef IDOCP_CONSTRAINTS_JOINT_VARIABLES_UPPER_LIMITS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_VARIABLES_UPPER_LIMITS_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "pdipm_test_func.hpp"


namespace idocp {
namespace pdipmtest {

class JointVariablesUpperLimits {
public:
  JointVariablesUpperLimits(const Robot& robot, const Eigen::VectorXd& xmax, 
                            const double barrier);

  JointVariablesUpperLimits();

  ~JointVariablesUpperLimits();

  // Use default copy constructor.
  JointVariablesUpperLimits(const JointVariablesUpperLimits&) = default;

  // Use default copy operator.
  JointVariablesUpperLimits& operator=(const JointVariablesUpperLimits&) 
      = default;

  // Use default move constructor.
  JointVariablesUpperLimits(JointVariablesUpperLimits&&) noexcept = default;

  // Use default move assign operator.
  JointVariablesUpperLimits& operator=(JointVariablesUpperLimits&&) noexcept
      = default;

  bool isFeasible(const Eigen::VectorXd& v);

  void setSlackAndDual(const double dtau, const Eigen::VectorXd& x);

  void condenseSlackAndDual(const double dtau, const Eigen::VectorXd& x, 
                            Eigen::MatrixXd& Cxx, Eigen::VectorXd& Cx);

  void computeSlackAndDualDirection(const double dtau,
                                    const Eigen::VectorXd& dx);

  double maxSlackStepSize(const double margin_rate);

  double maxDualStepSize(const double margin_rate);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  void augmentDualResidual(const double dtau, Eigen::VectorXd& Cx);

  double residualL1Nrom(const double dtau, const Eigen::VectorXd& x);

  double residualSquaredNrom(const double dtau, const Eigen::VectorXd& x);

private:
  int dimv_, dimc_, dim_passive_;
  bool has_floating_base_;
  double barrier_;
  Eigen::VectorXd xmax_, slack_, dual_, residual_, duality_, dslack_, ddual_;
};


inline JointVariablesUpperLimits::JointVariablesUpperLimits(
    const Robot& robot, const Eigen::VectorXd& xmax, const double barrier)
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


inline JointVariablesUpperLimits::JointVariablesUpperLimits()
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


inline bool JointVariablesUpperLimits::isFeasible(const Eigen::VectorXd& x) {
  assert(x.size() >= dimv_);
  for (int i=0; i<dimc_; ++i) {
    if (x.tail(dimc_).coeff(i) > xmax_.coeff(i)) {
      return false;
    }
  }
  return true;
}


inline void JointVariablesUpperLimits::setSlackAndDual(
    const double dtau, const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  slack_ = dtau * (xmax_-x.tail(dimc_));
  pdipmtestfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


inline void JointVariablesUpperLimits::condenseSlackAndDual(
    const double dtau, const Eigen::VectorXd& x, Eigen::MatrixXd& Cxx, 
    Eigen::VectorXd& Cx) {
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
  pdipmtestfunc::ComputeDualityResidual(barrier_, slack_, dual_, duality_);
  Cx.tail(dimc_).array() 
      += dtau * (dual_.array()*residual_.array()-duality_.array()) 
              / slack_.array();
}


inline void JointVariablesUpperLimits::computeSlackAndDualDirection(
    const double dtau, const Eigen::VectorXd& dx) {
  assert(dtau > 0);
  assert(dx.size() == dimv_);
  dslack_ = - dtau * dx.tail(dimc_) - residual_;
  pdipmtestfunc::ComputeDualDirection(dual_, slack_, dslack_, duality_, ddual_);
}


inline double JointVariablesUpperLimits::maxSlackStepSize(
    const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmtestfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


inline double JointVariablesUpperLimits::maxDualStepSize(
    const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmtestfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


inline void JointVariablesUpperLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


inline void JointVariablesUpperLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


inline double JointVariablesUpperLimits::costSlackBarrier() {
  return pdipmtestfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


inline double JointVariablesUpperLimits::costSlackBarrier(
    const double step_size) {
  return pdipmtestfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


inline void JointVariablesUpperLimits::augmentDualResidual(
    const double dtau, Eigen::VectorXd& Cx) {
  assert(dtau > 0);
  assert(Cx.size() == dimv_);
  Cx.tail(dimc_).noalias() += dtau * dual_;
}


inline double JointVariablesUpperLimits::residualL1Nrom(
    const double dtau, const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (x.tail(dimc_)-xmax_) + slack_;
  return residual_.lpNorm<1>();
}


inline double JointVariablesUpperLimits::residualSquaredNrom(
    const double dtau, const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (x.tail(dimc_)-xmax_) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  pdipmtestfunc::ComputeDualityResidual(barrier_, slack_, dual_, residual_);
  error += residual_.squaredNorm();
  return error;
}

} // namespace pdipmtest
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_JOINT_VARIABLES_UPPER_LIMITS_HPP_