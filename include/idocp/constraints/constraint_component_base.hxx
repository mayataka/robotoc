#ifndef IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_
#define IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_

#include "idocp/constraints/pdipm.hpp"
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

inline double ConstraintComponentBase::l1NormPrimalResidual(
    const ConstraintComponentData& data) const {
  return data.residual.lpNorm<1>();
}


inline double ConstraintComponentBase::squaredNormPrimalAndDualResidual(
    const ConstraintComponentData& data) const {
  return (data.residual.squaredNorm() + data.duality.squaredNorm());
}


inline double ConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundarySlack(fraction_to_boundary_rate_, data);
}


inline double ConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::FractionToBoundaryDual(fraction_to_boundary_rate_, data);
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
  return pdipm::CostBarrier(barrier_, data.slack);
}


inline double ConstraintComponentBase::costSlackBarrier(
    const ConstraintComponentData& data, const double step_size) const {
  return pdipm::CostBarrier(barrier_, data.slack+step_size*data.dslack);
}


inline void ConstraintComponentBase::setSlackAndDualPositive(
    ConstraintComponentData& data) const {
  pdipm::SetSlackAndDualPositive(barrier_, data);
}


inline void ConstraintComponentBase::computeDuality(
    ConstraintComponentData& data) const {
  pdipm::ComputeDuality(barrier_, data);
}


inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) const {
  pdipm::ComputeDualDirection(data);
}


inline double ConstraintComponentBase::computeDuality(const double slack, 
                                                      const double dual) const {
  return (slack * dual - barrier_); 
}


inline double ConstraintComponentBase::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double duality) const {
  return (- (dual * dslack + duality) / slack);
}

} // namespace idocp

#endif // IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_