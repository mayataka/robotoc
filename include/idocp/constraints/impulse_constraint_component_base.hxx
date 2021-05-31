#ifndef IDOCP_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_ 
#define IDOCP_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_

#include "idocp/constraints/impulse_constraint_component_base.hpp"
#include "idocp/constraints/pdipm.hpp"

#include <cassert>


namespace idocp {

inline double ImpulseConstraintComponentBase::l1NormPrimalResidual(
    const ConstraintComponentData& data) const {
  return data.residual.lpNorm<1>();
}


inline double ImpulseConstraintComponentBase::squaredNormPrimalAndDualResidual(
    const ConstraintComponentData& data) const {
  return (data.residual.squaredNorm() + data.duality.squaredNorm());
}


inline double ImpulseConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundarySlack(fraction_to_boundary_rate_, data);
}


inline double ImpulseConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundaryDual(fraction_to_boundary_rate_, data);
}


inline void ImpulseConstraintComponentBase::updateSlack(
    ConstraintComponentData& data, const double step_size) const {
  assert(step_size > 0);
  data.slack.noalias() += step_size * data.dslack;
}


inline void ImpulseConstraintComponentBase::updateDual(
    ConstraintComponentData& data, const double step_size) const {
  assert(step_size > 0);
  data.dual.noalias() += step_size * data.ddual;
}


inline double ImpulseConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data) const {
  return pdipm::CostBarrier(barrier_, data.slack);
}


inline double ImpulseConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data, const double step_size) const {
  return pdipm::CostBarrier(barrier_, data.slack+step_size*data.dslack);
}


inline double ImpulseConstraintComponentBase::barrier() const {
  return barrier_;
}


inline double ImpulseConstraintComponentBase::fractionToBoundaryRate() const {
  return fraction_to_boundary_rate_;
}


inline void ImpulseConstraintComponentBase::setBarrier(const double barrier) {
  assert(barrier > 0);
  barrier_ = barrier;
}


inline void ImpulseConstraintComponentBase::setFractionToBoundaryRate(
    const double fraction_to_boundary_rate) {
  assert(fraction_to_boundary_rate > 0);
  fraction_to_boundary_rate_ = fraction_to_boundary_rate;
}


inline void ImpulseConstraintComponentBase::setSlackAndDualPositive(
    ConstraintComponentData& data) const {
  pdipm::SetSlackAndDualPositive(barrier_, data);
}


inline void ImpulseConstraintComponentBase::computeDuality(
    ConstraintComponentData& data) const {
  pdipm::ComputeDuality(barrier_, data);
}


inline double ImpulseConstraintComponentBase::computeDuality(
    const double slack, const double dual) const {
  return (slack * dual - barrier_); 
}


inline void ImpulseConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) {
  pdipm::ComputeDualDirection(data);
}


inline double ImpulseConstraintComponentBase::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double duality) {
  return (- (dual * dslack + duality) / slack);
}


} // namespace idocp

#endif // IDOCP_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_ 