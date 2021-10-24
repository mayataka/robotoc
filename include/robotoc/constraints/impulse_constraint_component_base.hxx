#ifndef ROBOTOC_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_ 
#define ROBOTOC_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_

#include "robotoc/constraints/impulse_constraint_component_base.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

inline ImpulseConstraintComponentBase::ImpulseConstraintComponentBase(
    const double barrier, const double fraction_to_boundary_rule) 
  : barrier_(barrier),
    fraction_to_boundary_rule_(fraction_to_boundary_rule) {
    try {
    if (barrier <= 0) {
      throw std::out_of_range(
          "invalid argment: barrirer must be positive");
    }
    if (fraction_to_boundary_rule <= 0) {
      throw std::out_of_range(
          "invalid argment: fraction_to_boundary_rule must be positive");
    }
    if (fraction_to_boundary_rule >= 1) {
      throw std::out_of_range(
          "invalid argment: fraction_to_boundary_rule must be less than 1");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ImpulseConstraintComponentBase::ImpulseConstraintComponentBase() 
  : barrier_(0),
    fraction_to_boundary_rule_(0) {
}


inline double ImpulseConstraintComponentBase::maxSlackStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundarySlack(fraction_to_boundary_rule_, data);
}


inline double ImpulseConstraintComponentBase::maxDualStepSize(
    const ConstraintComponentData& data) const {
  return pdipm::fractionToBoundaryDual(fraction_to_boundary_rule_, data);
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
  pdipm::setSlackAndDualPositive(barrier_, data);
}


inline void ImpulseConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data) const {
  pdipm::computeComplementarySlackness(barrier_, data);
}


inline void ImpulseConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start, const int size) const {
  pdipm::computeComplementarySlackness(barrier_, data, start, size);
}


template <int Size>
inline void ImpulseConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start) const {
  pdipm::computeComplementarySlackness<Size>(barrier_, data, start);
}


inline double ImpulseConstraintComponentBase::computeComplementarySlackness(
    const double slack, const double dual) const {
  return pdipm::computeComplementarySlackness(barrier_, slack, dual);
}


inline void ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data) {
  pdipm::computeCondensingCoeffcient(data);
}


inline void ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::computeCondensingCoeffcient(data, start, size);
}


template <int Size>
inline void ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start) {
  pdipm::computeCondensingCoeffcient<Size>(data, start);
}


inline double ImpulseConstraintComponentBase::computeCondensingCoeffcient(
    const double slack, const double dual, const double residual, 
    const double cmpl) {
  return pdipm::computeCondensingCoeffcient(slack, dual, residual, cmpl);
}


inline void ImpulseConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) {
  pdipm::computeDualDirection(data);
}


inline void ImpulseConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::computeDualDirection(data, start, size);
}


template <int Size>
inline void ImpulseConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start) {
  pdipm::computeDualDirection<Size>(data, start);
}


inline double ImpulseConstraintComponentBase::computeDualDirection(
    const double slack, const double dual, const double dslack, 
    const double cmpl) {
  return (- (dual * dslack + cmpl) / slack);
}


template <typename VectorType>
inline double ImpulseConstraintComponentBase::logBarrier(
    const Eigen::MatrixBase<VectorType>& slack) const {
  return pdipm::logBarrier(barrier_, slack);
}

} // namespace robotoc

#endif // ROBOTOC_IMPULSE_CONSTRAINT_COMPONENT_BASE_HXX_ 