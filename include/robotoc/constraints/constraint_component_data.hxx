#ifndef ROBOTOC_CONSTRAINT_COMPONENT_DATA_HXX_
#define ROBOTOC_CONSTRAINT_COMPONENT_DATA_HXX_

#include "robotoc/constraints/constraint_component_data.hpp"

#include <cmath>
#include <stdexcept>
#include <iostream>


namespace robotoc {

inline ConstraintComponentData::ConstraintComponentData(const int dimc,
                                                        const double barrier_param)
  : slack(Eigen::VectorXd::Constant(dimc, std::sqrt(barrier_param))),
    dual(Eigen::VectorXd::Constant(dimc, std::sqrt(barrier_param))),
    residual(Eigen::VectorXd::Zero(dimc)),
    cmpl(Eigen::VectorXd::Zero(dimc)),
    dslack(Eigen::VectorXd::Zero(dimc)),
    ddual(Eigen::VectorXd::Zero(dimc)),
    cond(Eigen::VectorXd::Zero(dimc)),
    log_barrier(0),
    r(),
    J(),
    dimc_(dimc) {
  if (dimc <= 0) {
    throw std::out_of_range(
        "[ConstraintComponentData] invalid argment: 'dimc' must be positive!");
  }
  if (barrier_param <= 0) {
    throw std::out_of_range(
        "[ConstraintComponentData] invalid argment: 'barrier_param' must be positive!");
  }
}


inline ConstraintComponentData::ConstraintComponentData()
  : slack(),
    dual(),
    residual(),
    cmpl(),
    dslack(),
    ddual(),
    cond(),
    log_barrier(0),
    r(),
    J(),
    dimc_(0) {
}


inline void ConstraintComponentData::copySlackAndDual(
    const ConstraintComponentData& other) {
  assert(dimc() == other.dimc());
  slack = other.slack;
  dual = other.dual;
}


inline double ConstraintComponentData::KKTError() const {
  return (residual.squaredNorm() + cmpl.squaredNorm());
}


inline void ConstraintComponentData::resize(const int dimc) {
  assert(dimc >= 0);
  slack.conservativeResize(dimc);
  dual.conservativeResize(dimc);
  residual.resize(dimc);
  cmpl.resize(dimc);
  dslack.resize(dimc);
  ddual.resize(dimc);
  cond.resize(dimc);
}


inline int ConstraintComponentData::dimc() const {
  return dimc_;
}


inline bool ConstraintComponentData::checkDimensionalConsistency() const {
  if (slack.size() != dimc_) {
    return false;
  }
  if (dual.size() != dimc_) {
    return false;
  }
  if (residual.size() != dimc_) {
    return false;
  }
  if (cmpl.size() != dimc_) {
    return false;
  }
  if (dslack.size() != dimc_) {
    return false;
  }
  if (ddual.size() != dimc_) {
    return false;
  }
  if (cond.size() != dimc_) {
    return false;
  }
  return true;
}


inline bool ConstraintComponentData::isApprox(
    const ConstraintComponentData& other) const {
  if (!slack.isApprox(other.slack)) {
    return false;
  }
  if (!dual.isApprox(other.dual)) {
    return false;
  }
  if (!residual.isApprox(other.residual)) {
    return false;
  }
  if (!cmpl.isApprox(other.cmpl)) {
    return false;
  }
  if (!dslack.isApprox(other.dslack)) {
    return false;
  }
  if (!ddual.isApprox(other.ddual)) {
    return false;
  }
  if (!cond.isApprox(other.cond)) {
    return false;
  }
  Eigen::VectorXd lb(1), other_lb(1);
  lb << log_barrier;
  other_lb << other.log_barrier;
  if (!lb.isApprox(other_lb)) {
    return false;
  }
  return true;
}

} // namespace robotoc

#endif // ROBOTOC_CONSTRAINT_COMPONENT_DATA_HXX_