#ifndef ROBOTOC_DWELL_TIME_LOWER_BOUND_HXX_
#define ROBOTOC_DWELL_TIME_LOWER_BOUND_HXX_

#include "robotoc/hybrid/dwell_time_lower_bound.hpp"
#include "robotoc/constraints/pdipm.hpp"

#include <cmath>
#include <cassert>
#include <stdexcept>
#include <iostream>


namespace robotoc {

inline DwellTimeLowerBound::DwellTimeLowerBound(
    const double barrier_param, const double _fraction_to_boundary_rule) 
  : barrier_(barrier_param),
    fraction_to_boundary_rule_(_fraction_to_boundary_rule),
    slack_(std::sqrt(barrier_param)), 
    dual_(std::sqrt(barrier_param)), 
    residual_(0), 
    cmpl_(0), 
    dslack_(0), 
    ddual_(0), 
    log_barrier_(0) {
  if (barrier_param <= 0) {
    throw std::out_of_range(
        "[DwellTimeLowerBound] invalid argment: 'barrier_param' must be positive!");
  }
  if (_fraction_to_boundary_rule <= 0) {
    throw std::out_of_range(
        "[DwellTimeLowerBound] invalid argment: 'fraction_to_boundary_rule' must be positive!");
  }
  if (_fraction_to_boundary_rule >= 1) {
    throw std::out_of_range(
        "[DwellTimeLowerBound] invalid argment: 'fraction_to_boundary_rule' must be less than 1!");
  }
}


inline DwellTimeLowerBound::DwellTimeLowerBound() 
  : barrier_(0),
    fraction_to_boundary_rule_(0),
    slack_(0), 
    dual_(0), 
    residual_(0), 
    cmpl_(0), 
    dslack_(0), 
    ddual_(0), 
    log_barrier_(0) {
}


inline DwellTimeLowerBound::~DwellTimeLowerBound() {
}


inline void DwellTimeLowerBound::setSlack(const double min_dt, const double ts1, 
                                          const double ts2) {
  assert(min_dt >= 0.0);
  slack_ = ts2 - ts1 - min_dt;
  if (slack_ <= 0) {
    slack_ = std::sqrt(barrier_);
  }
  dual_ = barrier_ / slack_;
}


inline void DwellTimeLowerBound::evalConstraint(const double min_dt, 
                                                const double ts1, 
                                                const double ts2) {
  assert(min_dt >= 0.0);
  residual_ = ts1 + min_dt - ts2 + slack_;
  cmpl_  = pdipm::computeComplementarySlackness(barrier_, slack_, dual_);
  log_barrier_ = - barrier_ * std::log(slack_);
}


inline void DwellTimeLowerBound::evalDerivatives_lub(
    SplitKKTResidual& kkt_residual1, SplitKKTResidual& kkt_residual2) const {
  kkt_residual1.h -= dual_;
  kkt_residual2.h += dual_;
}


inline void DwellTimeLowerBound::evalDerivatives_ub(
    SplitKKTResidual& kkt_residual1) const {
  kkt_residual1.h -= dual_;
}


inline void DwellTimeLowerBound::evalDerivatives_lb(
    SplitKKTResidual& kkt_residual2) const {
  kkt_residual2.h += dual_;
}


inline void DwellTimeLowerBound::condenseSlackAndDual_lub(
    SplitKKTMatrix& kkt_matrix1, SplitKKTResidual& kkt_residual1,
    SplitKKTMatrix& kkt_matrix2, SplitKKTResidual& kkt_residual2) const {
  const double cond = 
      pdipm::computeCondensingCoeffcient(slack_, dual_, residual_, cmpl_);
  kkt_matrix1.Qtt += dual_ / slack_;
  kkt_residual1.h -= cond;
  kkt_matrix2.Qtt += dual_ / slack_;
  kkt_residual2.h += cond;
}


inline void DwellTimeLowerBound::condenseSlackAndDual_ub(
    SplitKKTMatrix& kkt_matrix1, SplitKKTResidual& kkt_residual1) const {
  const double cond = 
      pdipm::computeCondensingCoeffcient(slack_, dual_, residual_, cmpl_);
  kkt_matrix1.Qtt += dual_ / slack_;
  kkt_residual1.h -= cond;
}


inline void DwellTimeLowerBound::condenseSlackAndDual_lb(
    SplitKKTMatrix& kkt_matrix2, SplitKKTResidual& kkt_residual2) const {
  const double cond = 
      pdipm::computeCondensingCoeffcient(slack_, dual_, residual_, cmpl_);
  kkt_matrix2.Qtt += dual_ / slack_;
  kkt_residual2.h += cond;
}


inline void DwellTimeLowerBound::expandSlackAndDual_lub(const double dts1, 
                                                        const double dts2) {
  dslack_ = - dts1 + dts2 - residual_;
  ddual_  = - (dual_*dslack_+cmpl_) / slack_;
}


inline void DwellTimeLowerBound::expandSlackAndDual_ub(const double dts1) {
  dslack_ = - dts1 - residual_;
  ddual_  = - (dual_*dslack_+cmpl_) / slack_;
}


inline void DwellTimeLowerBound::expandSlackAndDual_lb(const double dts2) {
  dslack_ = dts2 - residual_;
  ddual_  = - (dual_*dslack_+cmpl_) / slack_;
}


inline double DwellTimeLowerBound::maxSlackStepSize() const {
  return pdipm::fractionToBoundary(fraction_to_boundary_rule_, slack_, dslack_);
}


inline double DwellTimeLowerBound::maxDualStepSize() const {
  return pdipm::fractionToBoundary(fraction_to_boundary_rule_, dual_, ddual_);
}


inline void DwellTimeLowerBound::updateSlack(const double step_size) {
  assert(step_size > 0.0);
  assert(step_size <= 1.0);
  slack_ += step_size * dslack_;
}


inline void DwellTimeLowerBound::updateDual(const double step_size) {
  assert(step_size > 0.0);
  assert(step_size <= 1.0);
  dual_ += step_size * ddual_;
}


inline double DwellTimeLowerBound::KKTError() const {
  return (residual_*residual_) + (cmpl_*cmpl_);
}


inline void DwellTimeLowerBound::setBarrierParam(const double barrier_param) {
  assert(barrier_param > 0.0);
  barrier_ = barrier_param;
}


inline void DwellTimeLowerBound::setFractionToBoundaryRule(
    const double _fraction_to_boundary_rule) {
  assert(_fraction_to_boundary_rule > 0.0);
  assert(_fraction_to_boundary_rule < 1.0);
  fraction_to_boundary_rule_ = _fraction_to_boundary_rule;
}


inline double DwellTimeLowerBound::getBarrierParam() const {
  return barrier_;
}


inline double DwellTimeLowerBound::getFractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}


inline void DwellTimeLowerBound::disp(std::ostream& os) const {
  os << "barrier_param = " << barrier_ << std::endl;
  os << "fraction_to_boundary_rule = " << fraction_to_boundary_rule_ << std::endl;
  os << "slack = " << slack_ << std::endl;
  os << "dual = " << dual_ << std::endl;
  os << "residual_ = " << residual_ << std::endl;
  os << "cmpl_ = " << cmpl_ << std::flush;
}


inline std::ostream& operator<<(std::ostream& os, 
                                const DwellTimeLowerBound& dtlb) {
  os << "DwellTimeLowerBound :" << std::endl;
  dtlb.disp(os);
  return os;
}

} // namespace robotoc

#endif // ROBOTOC_DWELL_TIME_LOWER_BOUND_HXX_ 