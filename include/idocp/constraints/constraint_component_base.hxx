#ifndef IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_
#define IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_

#include "idocp/constraints/pdipm_func.hpp"
#include <assert.h>


namespace idocp {

inline void ConstraintComponentBase::setBarrier(const double barrier) {
  assert(barrier > 0);
  barrier_ = barrier;
}


inline void ConstraintComponentBase::setFractionToBoundaryRate(
    const double fraction_to_boundary_rate) {
  assert(fraction_to_boundary_rate > 0);
  fraction_to_boundary_rate_ = fraction_to_boundary_rate;
}


inline double ConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return fractionToBoundary(data.slack, data.dslack);
}


inline double ConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return fractionToBoundary(data.dual, data.ddual);
}


inline void ConstraintComponentBase::updateSlack(ConstraintComponentData& data, 
                                                 const double step_size) const {
  assert(step_size > 0);
  data.slack.noalias() += step_size * data.dslack;
}


inline void ConstraintComponentBase::updateDual(ConstraintComponentData& data, 
                                                const double step_size) const {
  assert(step_size > 0);
  data.dual.noalias() += step_size * data.ddual;
}


inline double ConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data) const {
  const double cost = - barrier_ * data.slack.array().log().sum();
  return cost;
}


inline double ConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data, const double step_size) const {
  const double cost 
      = - barrier_ * (data.slack+step_size*data.dslack).array().log().sum();
  return cost;
}


inline void ConstraintComponentBase::setSlackAndDualPositive(
    Eigen::VectorXd& slack, Eigen::VectorXd& dual) const {
  pdipmfunc::SetSlackAndDualPositive(barrier_, slack, dual);
}


inline void ConstraintComponentBase::computeDuality(
    const Eigen::VectorXd& slack, const Eigen::VectorXd& dual, 
    Eigen::VectorXd& duality) const {
  pdipmfunc::ComputeDuality(barrier_, slack, dual, duality);
}


inline void ConstraintComponentBase::computeDualDirection(
    const Eigen::VectorXd& slack, const Eigen::VectorXd& dual, 
    const Eigen::VectorXd& dslack, const Eigen::VectorXd& duality, 
    Eigen::VectorXd& ddual) const {
  pdipmfunc::ComputeDualDirection(slack, dual, dslack, duality, ddual);
}


inline double ConstraintComponentBase::fractionToBoundary(
    const Eigen::VectorXd& vec, const Eigen::VectorXd& dvec) const {
  return pdipmfunc::FractionToBoundary(dimc(), fraction_to_boundary_rate_, 
                                       vec, dvec);
}

} // namespace idocp

#endif // IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_