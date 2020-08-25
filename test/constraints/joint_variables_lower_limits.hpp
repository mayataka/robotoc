#ifndef IDOCP_CONSTRAINTS_JOINT_VARIABLES_LOWER_LIMITS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_VARIABLES_LOWER_LIMITS_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "pdipm_test_func.hpp"


namespace idocp {
namespace pdipmtest {

class JointVariablesLowerLimits {
public:
  JointVariablesLowerLimits(const Robot& robot, const Eigen::VectorXd& xmin, 
                            const double barrier);

  JointVariablesLowerLimits();

  ~JointVariablesLowerLimits();

  // Use default copy constructor.
  JointVariablesLowerLimits(const JointVariablesLowerLimits&) = default;

  // Use default copy operator.
  JointVariablesLowerLimits& operator=(const JointVariablesLowerLimits&) 
      = default;

  // Use default move constructor.
  JointVariablesLowerLimits(JointVariablesLowerLimits&&) noexcept = default;

  // Use default move assign operator.
  JointVariablesLowerLimits& operator=(JointVariablesLowerLimits&&) noexcept 
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
  Eigen::VectorXd xmin_, slack_, dual_, residual_, duality_, dslack_, ddual_;
};


inline JointVariablesLowerLimits::JointVariablesLowerLimits(
    const Robot& robot, const Eigen::VectorXd& xmin, const double barrier)
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


inline JointVariablesLowerLimits::JointVariablesLowerLimits()
  : dimv_(0),
    dimc_(0),
    dim_passive_(0),
    has_floating_base_(false),
    barrier_(0),
    xmin_(),
    slack_(),
    dual_(),
    residual_(),
    duality_(),
    dslack_(), 
    ddual_() {
}


inline JointVariablesLowerLimits::~JointVariablesLowerLimits() {
}


inline bool JointVariablesLowerLimits::isFeasible(const Eigen::VectorXd& x) {
  assert(x.size() >= dimv_);
  for (int i=0; i<dimc_; ++i) {
    if (x.tail(dimc_).coeff(i) < xmin_.coeff(i)) {
      return false;
    }
  }
  return true;
}


inline void JointVariablesLowerLimits::setSlackAndDual(
    const double dtau, const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  slack_ = dtau * (x.tail(dimc_)-xmin_);
  pdipmtestfunc::SetSlackAndDualPositive(dimc_, barrier_, slack_, dual_);
}


inline void JointVariablesLowerLimits::condenseSlackAndDual(
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
  residual_ = dtau * (xmin_-x.tail(dimc_)) + slack_;
  pdipmtestfunc::ComputeDualityResidual(barrier_, slack_, dual_, duality_);
  Cx.tail(dimc_).array() 
      -= dtau * (dual_.array()*residual_.array()-duality_.array()) 
              / slack_.array();
}


inline void JointVariablesLowerLimits::computeSlackAndDualDirection(
    const double dtau, const Eigen::VectorXd& dx) {
  assert(dtau > 0);
  assert(dx.size() == dimv_);
  dslack_ = dtau * dx.tail(dimc_) - residual_;
  pdipmtestfunc::ComputeDualDirection(dual_, slack_, dslack_, duality_, ddual_);
}


inline double JointVariablesLowerLimits::maxSlackStepSize(
    const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmtestfunc::FractionToBoundary(dimc_, margin_rate, slack_, dslack_);
}


inline double JointVariablesLowerLimits::maxDualStepSize(
    const double margin_rate) {
  assert(margin_rate > 0);
  return pdipmtestfunc::FractionToBoundary(dimc_, margin_rate, dual_, ddual_);
}


inline void JointVariablesLowerLimits::updateSlack(const double step_size) {
  assert(step_size > 0);
  slack_.noalias() += step_size * dslack_;
}


inline void JointVariablesLowerLimits::updateDual(const double step_size) {
  assert(step_size > 0);
  dual_.noalias() += step_size * ddual_;
}


inline double JointVariablesLowerLimits::costSlackBarrier() {
  return pdipmtestfunc::SlackBarrierCost(dimc_, barrier_, slack_);
}


inline double JointVariablesLowerLimits::costSlackBarrier(const double step_size) {
  return pdipmtestfunc::SlackBarrierCost(dimc_, barrier_, slack_+step_size*dslack_);
}


inline void JointVariablesLowerLimits::augmentDualResidual(
    const double dtau, Eigen::VectorXd& Cx) {
  assert(dtau > 0);
  assert(Cx.size() == dimv_);
  Cx.tail(dimc_).noalias() -= dtau * dual_;
}


inline double JointVariablesLowerLimits::residualL1Nrom(
    const double dtau, const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (xmin_-x.tail(dimc_)) + slack_;
  return residual_.lpNorm<1>();
}


inline double JointVariablesLowerLimits::residualSquaredNrom(
    const double dtau, const Eigen::VectorXd& x) {
  assert(dtau > 0);
  assert(x.size() >= dimv_);
  residual_ = dtau * (xmin_-x.tail(dimc_)) + slack_;
  double error = 0;
  error += residual_.squaredNorm();
  pdipmtestfunc::ComputeDualityResidual(barrier_, slack_, dual_, residual_);
  error += residual_.squaredNorm();
  return error;
}


} // namespace pdipm
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_JOINT_VARIABLES_LOWER_LIMITS_HPP_