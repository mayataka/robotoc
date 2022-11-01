#ifndef ROBOTOC_CONSTRAINT_COMPONENT_BASE_HXX_
#define ROBOTOC_CONSTRAINT_COMPONENT_BASE_HXX_

#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include <cassert>
#include <stdexcept>


namespace robotoc {

inline double ConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundarySlack(fraction_to_boundary_rule_, data);
}


inline double ConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundaryDual(fraction_to_boundary_rule_, data);
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


template <typename Derived>
inline std::shared_ptr<Derived> ConstraintComponentBase::as() const {
  Derived* derived_ptr = dynamic_cast<Derived*>(this);
  if (derived_ptr == nullptr) {
    throw std::runtime_error("[ConstraintComponentBase] runtime error: failed in down-casting!");
  }
  return std::shared_ptr<Derived>(derived_ptr);
}


inline void ConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data) const {
  pdipm::computeComplementarySlackness(barrier_, data);
}


inline void ConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start, const int size) const {
  pdipm::computeComplementarySlackness(barrier_, data, start, size);
}


template <int Size>
inline void ConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start) const {
  pdipm::computeComplementarySlackness<Size>(barrier_, data, start);
}


inline double ConstraintComponentBase::computeComplementarySlackness(
    const double slack, const double dual) const {
  return pdipm::computeComplementarySlackness(barrier_, slack, dual);
}


inline void ConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data) {
  pdipm::computeCondensingCoeffcient(data);
}


inline void ConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::computeCondensingCoeffcient(data, start, size);
}


template <int Size>
inline void ConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start) {
  pdipm::computeCondensingCoeffcient<Size>(data, start);
}


inline double ConstraintComponentBase::computeCondensingCoeffcient(
    const double slack, const double dual, const double residual, 
    const double cmpl) {
  return pdipm::computeCondensingCoeffcient(slack, dual, residual, cmpl);
}


inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) {
  pdipm::computeDualDirection(data);
}


inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::computeDualDirection(data, start, size);
}


template <int Size>
inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start) {
  pdipm::computeDualDirection<Size>(data, start);
}


inline double ConstraintComponentBase::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double cmpl) {
  return pdipm::computeDualDirection(slack, dual, dslack, cmpl);
}


template <typename VectorType>
inline double ConstraintComponentBase::logBarrier(
    const Eigen::MatrixBase<VectorType>& slack) const {
  return pdipm::logBarrier(barrier_, slack);
}

} // namespace robotoc

#endif // ROBOTOC_CONSTRAINT_COMPONENT_BASE_HXX_