#ifndef ROBOTOC_IMPACT_CONSTRAINT_COMPONENT_BASE_HXX_ 
#define ROBOTOC_IMPACT_CONSTRAINT_COMPONENT_BASE_HXX_

#include "robotoc/constraints/impact_constraint_component_base.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include <cassert>


namespace robotoc {

inline double ImpactConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundarySlack(fraction_to_boundary_rule_, data);
}


inline double ImpactConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundaryDual(fraction_to_boundary_rule_, data);
}


inline void ImpactConstraintComponentBase::updateSlack(
    ConstraintComponentData& data, const double step_size) {
  assert(step_size > 0);
  data.slack.noalias() += step_size * data.dslack;
}


inline void ImpactConstraintComponentBase::updateDual(
    ConstraintComponentData& data, const double step_size) {
  assert(step_size > 0);
  data.dual.noalias() += step_size * data.ddual;
}


template <typename Derived>
inline std::shared_ptr<Derived> ImpactConstraintComponentBase::as_shared_ptr() {
  auto ptr = shared_from_this();
  auto derived_ptr = std::dynamic_pointer_cast<Derived>(ptr);
  if (derived_ptr == nullptr) {
    throw std::runtime_error("[ImpactConstraintComponentBase] runtime error: failed in down-casting!");
  }
  return derived_ptr;
}


inline void ImpactConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data) const {
  pdipm::computeComplementarySlackness(barrier_, data);
}


inline void ImpactConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start, const int size) const {
  pdipm::computeComplementarySlackness(barrier_, data, start, size);
}


template <int Size>
inline void ImpactConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start) const {
  pdipm::computeComplementarySlackness<Size>(barrier_, data, start);
}


inline double ImpactConstraintComponentBase::computeComplementarySlackness(
    const double slack, const double dual) const {
  return pdipm::computeComplementarySlackness(barrier_, slack, dual);
}


inline void ImpactConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data) {
  pdipm::computeCondensingCoeffcient(data);
}


inline void ImpactConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::computeCondensingCoeffcient(data, start, size);
}


template <int Size>
inline void ImpactConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start) {
  pdipm::computeCondensingCoeffcient<Size>(data, start);
}


inline double ImpactConstraintComponentBase::computeCondensingCoeffcient(
    const double slack, const double dual, const double residual, 
    const double cmpl) {
  return pdipm::computeCondensingCoeffcient(slack, dual, residual, cmpl);
}


inline void ImpactConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) {
  pdipm::computeDualDirection(data);
}


inline void ImpactConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::computeDualDirection(data, start, size);
}


template <int Size>
inline void ImpactConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start) {
  pdipm::computeDualDirection<Size>(data, start);
}


inline double ImpactConstraintComponentBase::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double cmpl) {
  return (- (dual * dslack + cmpl) / slack);
}


template <typename VectorType>
inline double ImpactConstraintComponentBase::logBarrier(
    const Eigen::MatrixBase<VectorType>& slack) const {
  return pdipm::logBarrier(barrier_, slack);
}

} // namespace robotoc

#endif // ROBOTOC_IMPACT_CONSTRAINT_COMPONENT_BASE_HXX_ 