#include "robotoc/constraints/impact_constraint_component_base.hpp"

#include "robotoc/constraints/pdipm.hpp"

#include <stdexcept>
#include <iostream>


namespace robotoc {

ImpactConstraintComponentBase::ImpactConstraintComponentBase(
    const double barrier_param, const double fraction_to_boundary_rule) 
  : barrier_(barrier_param),
    fraction_to_boundary_rule_(fraction_to_boundary_rule) {
  if (barrier_param <= 0) {
    throw std::out_of_range(
        "[ImpactConstraintComponentBase] invalid argment: 'barrier_param' must be positive");
  }
  if (fraction_to_boundary_rule <= 0) {
    throw std::out_of_range(
        "[ImpactConstraintComponentBase] invalid argment: 'fraction_to_boundary_rule' must be positive");
  }
  if (fraction_to_boundary_rule >= 1) {
    throw std::out_of_range(
        "[ImpactConstraintComponentBase] invalid argment: 'fraction_to_boundary_rule' must be less than 1");
  }
}


double ImpactConstraintComponentBase::getBarrierParam() const {
  return barrier_;
}


double ImpactConstraintComponentBase::getFractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}


void ImpactConstraintComponentBase::setBarrierParam(const double barrier_param) {
  if (barrier_param <= 0) {
    throw std::out_of_range(
        "[ImpactConstraintComponentBase] invalid argment: 'barrier_param' must be positive");
  }
  barrier_ = barrier_param;
}


void ImpactConstraintComponentBase::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  if (fraction_to_boundary_rule <= 0) {
    throw std::out_of_range(
        "[ImpactConstraintComponentBase] invalid argment: 'fraction_to_boundary_rule' must be positive");
  }
  if (fraction_to_boundary_rule >= 1) {
    throw std::out_of_range(
        "[ImpactConstraintComponentBase] invalid argment: 'fraction_to_boundary_rule' must be less than 1");
  }
  fraction_to_boundary_rule_ = fraction_to_boundary_rule;
}


void ImpactConstraintComponentBase::setSlackAndDualPositive(
    ConstraintComponentData& data) const {
  pdipm::setSlackAndDualPositive(barrier_, data);
}

} // namespace robotoc