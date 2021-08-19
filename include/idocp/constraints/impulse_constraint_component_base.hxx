#ifndef IDOCP_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_ 
#define IDOCP_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_

#include "idocp/constraints/impulse_constraint_component_base.hpp"
#include "idocp/constraints/pdipm.hpp"

#include <cassert>


namespace idocp {

inline double ImpulseConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundarySlack(fraction_to_boundary_rule_, data);
}


inline double ImpulseConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundaryDual(fraction_to_boundary_rule_, data);
}


inline void ImpulseConstraintComponentBase::updateSlack(
    ConstraintComponentData& data, const double step_size) {
  assert(step_size > 0);
  data.slack.noalias() += step_size * data.dslack;
}


inline void ImpulseConstraintComponentBase::updateDual(
    ConstraintComponentData& data, const double step_size) {
  assert(step_size > 0);
  data.dual.noalias() += step_size * data.ddual;
}


inline double ImpulseConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data) const {
  return pdipm::LogBarrier(barrier_, data.slack);
}


inline double ImpulseConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data, const double step_size) const {
  return pdipm::LogBarrier(barrier_, data.slack+step_size*data.dslack);
}


inline double ImpulseConstraintComponentBase::barrierParameter() const {
  return barrier_;
}


inline double ImpulseConstraintComponentBase::fractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}


inline void ImpulseConstraintComponentBase::setBarrier(const double barrier) {
  assert(barrier > 0);
  barrier_ = barrier;
}


inline void ImpulseConstraintComponentBase::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  assert(fraction_to_boundary_rule > 0);
  fraction_to_boundary_rule_ = fraction_to_boundary_rule;
}


inline void ImpulseConstraintComponentBase::setSlackAndDualPositive(
    ConstraintComponentData& data) const {
  pdipm::SetSlackAndDualPositive(barrier_, data);
}


inline void ImpulseConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data) const {
  pdipm::ComputeComplementarySlackness(barrier_, data);
}


inline void ImpulseConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start, const int size) const {
  pdipm::ComputeComplementarySlackness(barrier_, data, start, size);
}


template <int Size>
inline void ImpulseConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start) const {
  pdipm::ComputeComplementarySlackness<Size>(barrier_, data, start);
}


inline double ImpulseConstraintComponentBase::computeComplementarySlackness(
    const double slack, const double dual) const {
  return pdipm::ComputeComplementarySlackness(barrier_, slack, dual);
}


inline void ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data) {
  pdipm::ComputeCondensingCoeffcient(data);
}


inline void ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::ComputeCondensingCoeffcient(data, start, size);
}


template <int Size>
inline void ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start) {
  pdipm::ComputeCondensingCoeffcient<Size>(data, start);
}


inline double ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    const double slack, const double dual, const double residual, 
    const double cmpl) {
  return pdipm::ComputeCondensingCoeffcient(slack, dual, residual, cmpl);
}


inline void ImpulseConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) {
  pdipm::ComputeDualDirection(data);
}


inline void ImpulseConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::ComputeDualDirection(data, start, size);
}


template <int Size>
inline void ImpulseConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start) {
  pdipm::ComputeDualDirection<Size>(data, start);
}


inline double ImpulseConstraintComponentBase::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double cmpl) {
  return (- (dual * dslack + cmpl) / slack);
}


template <typename VectorType>
inline double ImpulseConstraintComponentBase::logBarrier(
    const Eigen::MatrixBase<VectorType>& slack) const {
  return pdipm::LogBarrier(barrier_, slack);
}

} // namespace idocp

#endif // IDOCP_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_ 