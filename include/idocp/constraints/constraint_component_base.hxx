#ifndef IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_
#define IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_

#include "idocp/constraints/pdipm.hpp"
#include <cassert>


namespace idocp {

inline double ConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundarySlack(fraction_to_boundary_rule_, data);
}


inline double ConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundaryDual(fraction_to_boundary_rule_, data);
}


inline void ConstraintComponentBase::updateSlack(ConstraintComponentData& data, 
                                                 const double step_size) {
  assert(step_size > 0);
  data.slack.noalias() += step_size * data.dslack;
}


inline void ConstraintComponentBase::updateDual(ConstraintComponentData& data, 
                                                const double step_size) {
  assert(step_size > 0);
  data.dual.noalias() += step_size * data.ddual;
}


inline double ConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data) const {
  return pdipm::LogBarrier(barrier_, data.slack);
}


inline double ConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data, const double step_size) const {
  return pdipm::LogBarrier(barrier_, data.slack+step_size*data.dslack);
}


inline double ConstraintComponentBase::barrierParameter() const {
  return barrier_;
}


inline double ConstraintComponentBase::fractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}


inline void ConstraintComponentBase::setBarrier(const double barrier) {
  assert(barrier > 0);
  barrier_ = barrier;
}


inline void ConstraintComponentBase::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  assert(fraction_to_boundary_rule > 0);
  fraction_to_boundary_rule_ = fraction_to_boundary_rule;
}


inline void ConstraintComponentBase::setSlackAndDualPositive(
    ConstraintComponentData& data) const {
  pdipm::SetSlackAndDualPositive(barrier_, data);
}


inline void ConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data) const {
  pdipm::ComputeComplementarySlackness(barrier_, data);
}


inline void ConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start, const int size) const {
  pdipm::ComputeComplementarySlackness(barrier_, data, start, size);
}


inline double ConstraintComponentBase::computeComplementarySlackness(
    const double slack, const double dual) const {
  return pdipm::ComputeComplementarySlackness(barrier_, slack, dual);
}


inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) {
  pdipm::ComputeDualDirection(data);
}


inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::ComputeDualDirection(data, start, size);
}


inline double ConstraintComponentBase::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double cmpl) {
  return pdipm::ComputeDualDirection(slack, dual, dslack, cmpl);
}


template <typename VectorType>
inline double ConstraintComponentBase::logBarrier(
    const Eigen::MatrixBase<VectorType>& slack) const {
  return pdipm::LogBarrier(barrier_, slack);
}

} // namespace idocp

#endif // IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_