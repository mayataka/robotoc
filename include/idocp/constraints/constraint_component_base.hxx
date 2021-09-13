#ifndef IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_
#define IDOCP_CONSTRAINT_COMPONENT_BASE_HXX_

#include "idocp/constraints/pdipm.hpp"
#include <cassert>
#include <stdexcept>
#include <iostream>


namespace idocp {

inline ConstraintComponentBase::ConstraintComponentBase(
    const double barrier, const double fraction_to_boundary_rule) 
  : barrier_(barrier),
    fraction_to_boundary_rule_(fraction_to_boundary_rule) {
  try {
    if (barrier <= 0) {
      throw std::out_of_range(
          "Invalid argment: barrirer must be positive!");
    }
    if (fraction_to_boundary_rule <= 0) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rule must be positive!");
    }
    if (fraction_to_boundary_rule >= 1) {
      throw std::out_of_range(
          "Invalid argment: fraction_to_boundary_rule must be less than 1!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ConstraintComponentBase::ConstraintComponentBase() 
  : barrier_(0),
    fraction_to_boundary_rule_(0) {
}


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


template <int Size>
inline void ConstraintComponentBase::computeComplementarySlackness(
    ConstraintComponentData& data, const int start) const {
  pdipm::ComputeComplementarySlackness<Size>(barrier_, data, start);
}


inline double ConstraintComponentBase::computeComplementarySlackness(
    const double slack, const double dual) const {
  return pdipm::ComputeComplementarySlackness(barrier_, slack, dual);
}


inline void ConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data) {
  pdipm::ComputeCondensingCoeffcient(data);
}


inline void ConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::ComputeCondensingCoeffcient(data, start, size);
}


template <int Size>
inline void ConstraintComponentBase::computeCondensingCoeffcient(
    ConstraintComponentData& data, const int start) {
  pdipm::ComputeCondensingCoeffcient<Size>(data, start);
}


inline double ConstraintComponentBase::computeCondensingCoeffcient(
    const double slack, const double dual, const double residual, 
    const double cmpl) {
  return pdipm::ComputeCondensingCoeffcient(slack, dual, residual, cmpl);
}


inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data) {
  pdipm::ComputeDualDirection(data);
}


inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start, const int size) {
  pdipm::ComputeDualDirection(data, start, size);
}


template <int Size>
inline void ConstraintComponentBase::computeDualDirection(
    ConstraintComponentData& data, const int start) {
  pdipm::ComputeDualDirection<Size>(data, start);
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