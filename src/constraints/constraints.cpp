#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/constraints_impl.hpp"

#include <stdexcept>
#include <cassert>


namespace robotoc {

Constraints::Constraints(const double barrier_param, 
                         const double fraction_to_boundary_rule) 
  : position_level_constraints_(),
    velocity_level_constraints_(),
    acceleration_level_constraints_(),
    impact_level_constraints_(),
    position_level_constraint_names_(), 
    velocity_level_constraint_names_(), 
    acceleration_level_constraint_names_(),
    impact_level_constraint_names_(),
    barrier_(barrier_param), 
    fraction_to_boundary_rule_(fraction_to_boundary_rule) {
  if (barrier_param <= 0) {
    throw std::out_of_range(
        "[Constraints] invalid argment: 'barrier_param' must be positive!");
  }
  if (fraction_to_boundary_rule <= 0) {
    throw std::out_of_range(
        "[Constraints] invalid argment: 'fraction_to_boundary_rule' must be positive!");
  }
  if (fraction_to_boundary_rule >= 1) {
    throw std::out_of_range(
        "[Constraints] invalid argment: 'fraction_to_boundary_rule' must be less than 1!");
  }
}


bool Constraints::exist(const std::string& name) const {
  return ((position_level_constraint_names_.find(name) != position_level_constraint_names_.end())
          || (velocity_level_constraint_names_.find(name) != velocity_level_constraint_names_.end())
          || (acceleration_level_constraint_names_.find(name) != acceleration_level_constraint_names_.end())
          || (impact_level_constraint_names_.find(name) != impact_level_constraint_names_.end()));
}


void Constraints::add(const std::string& name, 
                      ConstraintComponentBasePtr constraint) {
  if (exist(name)) {
    throw std::runtime_error("[Constraints] invalid argument: constraint component '" + name + "' already exists!");
  }
  constraint->setBarrierParam(barrier_);
  constraint->setFractionToBoundaryRule(fraction_to_boundary_rule_);
  if (constraint->kinematicsLevel() 
        == KinematicsLevel::PositionLevel) {
    position_level_constraint_names_.emplace(name, position_level_constraints_.size());
    position_level_constraints_.push_back(constraint);
  }
  else if (constraint->kinematicsLevel() 
              == KinematicsLevel::VelocityLevel) {
    velocity_level_constraint_names_.emplace(name, velocity_level_constraints_.size());
    velocity_level_constraints_.push_back(constraint);
  }
  else if (constraint->kinematicsLevel() 
              == KinematicsLevel::AccelerationLevel) {
    acceleration_level_constraint_names_.emplace(name, acceleration_level_constraints_.size());
    acceleration_level_constraints_.push_back(constraint);
  }
}


void Constraints::add(const std::string& name, 
                      ImpactConstraintComponentBasePtr constraint) {
  if (exist(name)) {
    throw std::runtime_error("[Constraints] invalid argument: constraint component '" + name + "' already exists!");
  }
  constraint->setBarrierParam(barrier_);
  constraint->setFractionToBoundaryRule(fraction_to_boundary_rule_);
  if (constraint->kinematicsLevel() 
              == KinematicsLevel::AccelerationLevel) {
    impact_level_constraint_names_.emplace(name, impact_level_constraints_.size());
    impact_level_constraints_.push_back(constraint);
  }
}


void Constraints::erase(const std::string& name) {
  if (!exist(name)) {
    throw std::runtime_error("[Constraints] invalid argument: constraint component '" + name + "' does not exist!");
  }
  if (position_level_constraint_names_.find(name) != position_level_constraint_names_.end()) {
    const int index = position_level_constraint_names_.at(name);
    position_level_constraint_names_.erase(name);
    position_level_constraints_.erase(position_level_constraints_.begin()+index);
  }
  else if (velocity_level_constraint_names_.find(name) != velocity_level_constraint_names_.end()) {
    const int index = velocity_level_constraint_names_.at(name);
    velocity_level_constraint_names_.erase(name);
    velocity_level_constraints_.erase(velocity_level_constraints_.begin()+index);
  }
  else if (acceleration_level_constraint_names_.find(name) != acceleration_level_constraint_names_.end()) {
    const int index = acceleration_level_constraint_names_.at(name);
    acceleration_level_constraint_names_.erase(name);
    acceleration_level_constraints_.erase(acceleration_level_constraints_.begin()+index);
  }
  else if (impact_level_constraint_names_.find(name) != impact_level_constraint_names_.end()) {
    const int index = impact_level_constraint_names_.at(name);
    impact_level_constraint_names_.erase(name);
    impact_level_constraints_.erase(impact_level_constraints_.begin()+index);
  }
}


std::shared_ptr<ConstraintComponentBase> 
Constraints::getConstraintComponent(const std::string& name) const {
  if (!exist(name)) {
    throw std::runtime_error("[Constraints] invalid argument: constraint component '" + name + "' does not exist!");
  }
  if (position_level_constraint_names_.find(name) != position_level_constraint_names_.end()) {
    const int index = position_level_constraint_names_.at(name);
    return position_level_constraints_.at(index);
  }
  else if (velocity_level_constraint_names_.find(name) != velocity_level_constraint_names_.end()) {
    const int index = velocity_level_constraint_names_.at(name);
    return velocity_level_constraints_.at(index);
  }
  else if (acceleration_level_constraint_names_.find(name) != acceleration_level_constraint_names_.end()) {
    const int index = acceleration_level_constraint_names_.at(name);
    return acceleration_level_constraints_.at(index);
  }
  else {
    throw std::runtime_error("[Constraints] invalid argument: constraint component '" + name + "' does not exist!");
  }
}


std::shared_ptr<ImpactConstraintComponentBase> 
Constraints::getImpactConstraintComponent(const std::string& name) const {
  if (!exist(name)) {
    throw std::runtime_error("[Constraints] invalid argument: constraint component '" + name + "' does not exist!");
  }
  if (impact_level_constraint_names_.find(name) != impact_level_constraint_names_.end()) {
    const int index = impact_level_constraint_names_.at(name);
    return impact_level_constraints_.at(index);
  }
  else {
    throw std::runtime_error("[Constraints] invalid argument: constraint component '" + name + "' does not exist!");
  }
}


void Constraints::clear() {
  position_level_constraints_.clear();
  velocity_level_constraints_.clear();
  acceleration_level_constraints_.clear();
  impact_level_constraints_.clear();
  position_level_constraint_names_.clear();
  velocity_level_constraint_names_.clear(); 
  acceleration_level_constraint_names_.clear();
  impact_level_constraint_names_.clear();
}


ConstraintsData Constraints::createConstraintsData(const Robot& robot, 
                                                   const int time_stage) const {
  ConstraintsData data(time_stage);
  constraintsimpl::createConstraintsData(position_level_constraints_, 
                                         data.position_level_data);
  constraintsimpl::createConstraintsData(velocity_level_constraints_, 
                                         data.velocity_level_data);
  constraintsimpl::createConstraintsData(acceleration_level_constraints_, 
                                         data.acceleration_level_data);
  constraintsimpl::createConstraintsData(impact_level_constraints_, 
                                         data.impact_level_data);
  return data;
}


bool Constraints::isFeasible(Robot& robot, const ContactStatus& contact_status, 
                             ConstraintsData& data, 
                             const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    if (!constraintsimpl::isFeasible(position_level_constraints_, robot, 
                                     contact_status, 
                                     data.position_level_data, s)) {
      return false;
    }
  }
  if (data.isVelocityLevelValid()) {
    if (!constraintsimpl::isFeasible(velocity_level_constraints_, robot, 
                                     contact_status, 
                                     data.velocity_level_data, s)) {
      return false;
    }
  }
  if (data.isAccelerationLevelValid()) {
    if (!constraintsimpl::isFeasible(acceleration_level_constraints_, robot, 
                                     contact_status, 
                                     data.acceleration_level_data, s)) {
      return false;
    }
  }
  return true;
}


bool Constraints::isFeasible(Robot& robot, const ImpactStatus& impact_status, 
                             ConstraintsData& data, 
                             const SplitSolution& s) const {
  if (data.isImpactLevelValid()) {
    if (!constraintsimpl::isFeasible(impact_level_constraints_, robot, 
                                     impact_status, 
                                     data.impact_level_data, s)) {
      return false;
    }
  }
  return true;
}


void Constraints::setSlackAndDual(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  ConstraintsData& data, 
                                  const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::setSlackAndDual(position_level_constraints_, robot, 
                                     contact_status, 
                                     data.position_level_data, s);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::setSlackAndDual(velocity_level_constraints_, robot, 
                                     contact_status, 
                                     data.velocity_level_data, s);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::setSlackAndDual(acceleration_level_constraints_, robot,
                                     contact_status, 
                                     data.acceleration_level_data, s);
  }
}


void Constraints::setSlackAndDual(Robot& robot, 
                                  const ImpactStatus& impact_status,
                                  ConstraintsData& data, 
                                  const SplitSolution& s) const {
  if (data.isImpactLevelValid()) {
    constraintsimpl::setSlackAndDual(impact_level_constraints_, robot,
                                     impact_status, 
                                     data.impact_level_data, s);
  }
}


void Constraints::evalConstraint(Robot& robot, 
                                 const ContactStatus& contact_status, 
                                 ConstraintsData& data, 
                                 const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::evalConstraint(position_level_constraints_, robot, 
                                    contact_status, data.position_level_data, s);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::evalConstraint(velocity_level_constraints_, robot, 
                                    contact_status, data.velocity_level_data, s);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::evalConstraint(acceleration_level_constraints_, robot, 
                                    contact_status, data.acceleration_level_data, s);
  }
}


void Constraints::evalConstraint(Robot& robot, 
                                 const ImpactStatus& impact_status, 
                                 ConstraintsData& data, 
                                 const SplitSolution& s) const {
  if (data.isImpactLevelValid()) {
    constraintsimpl::evalConstraint(impact_level_constraints_, robot, 
                                    impact_status, data.impact_level_data, s);
  }
}


void Constraints::linearizeConstraints(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       ConstraintsData& data, 
                                       const SplitSolution& s, 
                                       SplitKKTResidual& kkt_residual) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::linearizeConstraints(position_level_constraints_, robot, 
                                          contact_status, 
                                          data.position_level_data, s, 
                                          kkt_residual);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::linearizeConstraints(velocity_level_constraints_, robot, 
                                          contact_status, 
                                          data.velocity_level_data, s, 
                                          kkt_residual);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::linearizeConstraints(acceleration_level_constraints_, robot, 
                                          contact_status, 
                                          data.acceleration_level_data, 
                                          s, kkt_residual);
  }
}


void Constraints::linearizeConstraints(Robot& robot, 
                                       const ImpactStatus& impact_status, 
                                       ConstraintsData& data, 
                                       const SplitSolution& s, 
                                       SplitKKTResidual& kkt_residual) const {
  if (data.isImpactLevelValid()) {
    constraintsimpl::linearizeConstraints(impact_level_constraints_, robot, 
                                          impact_status, 
                                          data.impact_level_data, s, kkt_residual);
  }
}


void Constraints::condenseSlackAndDual(const ContactStatus& contact_status, 
                                       ConstraintsData& data, 
                                       SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::condenseSlackAndDual(position_level_constraints_, 
                                          contact_status, 
                                          data.position_level_data, 
                                          kkt_matrix, kkt_residual);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::condenseSlackAndDual(velocity_level_constraints_, 
                                          contact_status, 
                                          data.velocity_level_data, 
                                          kkt_matrix, kkt_residual);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::condenseSlackAndDual(acceleration_level_constraints_, 
                                          contact_status, 
                                          data.acceleration_level_data, 
                                          kkt_matrix, kkt_residual);
  }
}


void Constraints::condenseSlackAndDual(const ImpactStatus& impact_status, 
                                       ConstraintsData& data, 
                                       SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual) const {
  if (data.isImpactLevelValid()) {
    constraintsimpl::condenseSlackAndDual(impact_level_constraints_, 
                                          impact_status, 
                                          data.impact_level_data, 
                                          kkt_matrix, kkt_residual);
  }
}


void Constraints::expandSlackAndDual(const ContactStatus& contact_status,
                                     ConstraintsData& data, 
                                     const SplitDirection& d) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::expandSlackAndDual(position_level_constraints_, 
                                        contact_status,
                                        data.position_level_data, d);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::expandSlackAndDual(velocity_level_constraints_, 
                                        contact_status,
                                        data.velocity_level_data, d);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::expandSlackAndDual(acceleration_level_constraints_, 
                                        contact_status,
                                        data.acceleration_level_data, d);
  }
}


void Constraints::expandSlackAndDual(const ImpactStatus& impact_status, 
                                     ConstraintsData& data, 
                                     const SplitDirection& d) const {
  if (data.isImpactLevelValid()) {
    constraintsimpl::expandSlackAndDual(impact_level_constraints_, 
                                        impact_status, 
                                        data.impact_level_data, d);
  }
}


double Constraints::maxSlackStepSize(const ConstraintsData& data) const {
  double min_step_size = 1;
  if (data.isPositionLevelValid()) {
    const double step_size = constraintsimpl::maxSlackStepSize(
        position_level_constraints_, data.position_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isVelocityLevelValid()) {
    const double step_size = constraintsimpl::maxSlackStepSize(
        velocity_level_constraints_, data.velocity_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isAccelerationLevelValid()) {
    const double step_size = constraintsimpl::maxSlackStepSize(
        acceleration_level_constraints_, data.acceleration_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isImpactLevelValid()) {
    const double step_size = constraintsimpl::maxSlackStepSize(
        impact_level_constraints_, data.impact_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


double Constraints::maxDualStepSize(const ConstraintsData& data) const {
  double min_step_size = 1;
  if (data.isPositionLevelValid()) {
    const double step_size = constraintsimpl::maxDualStepSize(
        position_level_constraints_, data.position_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isVelocityLevelValid()) {
    const double step_size = constraintsimpl::maxDualStepSize(
        velocity_level_constraints_, data.velocity_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isAccelerationLevelValid()) {
    const double step_size = constraintsimpl::maxDualStepSize(
        acceleration_level_constraints_, data.acceleration_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isImpactLevelValid()) {
    const double step_size = constraintsimpl::maxDualStepSize(
        impact_level_constraints_, data.impact_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


void Constraints::updateSlack(ConstraintsData& data, const double step_size) {
  assert(step_size >= 0);
  assert(step_size <= 1);
  if (data.isPositionLevelValid()) {
    constraintsimpl::updateSlack(data.position_level_data, step_size);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::updateSlack(data.velocity_level_data, step_size);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::updateSlack(data.acceleration_level_data, step_size);
  }
  if (data.isImpactLevelValid()) {
    constraintsimpl::updateSlack(data.impact_level_data, step_size);
  }
}


void Constraints::updateDual(ConstraintsData& data, const double step_size) {
  assert(step_size >= 0);
  assert(step_size <= 1);
  if (data.isPositionLevelValid()) {
    constraintsimpl::updateDual(data.position_level_data, step_size);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::updateDual(data.velocity_level_data, step_size);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::updateDual(data.acceleration_level_data, step_size);
  }
  if (data.isImpactLevelValid()) {
    constraintsimpl::updateDual(data.impact_level_data, step_size);
  }
}


void Constraints::setBarrierParam(const double barrier_param) {
  if (barrier_param <= 0) {
    throw std::out_of_range(
        "[Constraints] invalid argment: 'barrier_param' must be positive");
  }
  constraintsimpl::setBarrierParam(position_level_constraints_, barrier_param);
  constraintsimpl::setBarrierParam(velocity_level_constraints_, barrier_param);
  constraintsimpl::setBarrierParam(acceleration_level_constraints_, barrier_param);
  constraintsimpl::setBarrierParam(impact_level_constraints_, barrier_param);
  barrier_ = barrier_param;
}


void Constraints::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  if (fraction_to_boundary_rule <= 0) {
    throw std::out_of_range(
        "[Constraint] invalid argment: 'fraction_to_boundary_rule' must be positive");
  }
  if (fraction_to_boundary_rule >= 1) {
    throw std::out_of_range(
        "[Constraint] invalid argment: 'fraction_to_boundary_rule' must be less than 1");
  }
  constraintsimpl::setFractionToBoundaryRule(position_level_constraints_, 
                                             fraction_to_boundary_rule);
  constraintsimpl::setFractionToBoundaryRule(velocity_level_constraints_, 
                                             fraction_to_boundary_rule);
  constraintsimpl::setFractionToBoundaryRule(acceleration_level_constraints_, 
                                             fraction_to_boundary_rule);
  constraintsimpl::setFractionToBoundaryRule(impact_level_constraints_, 
                                             fraction_to_boundary_rule);
  fraction_to_boundary_rule_ = fraction_to_boundary_rule;
}


double Constraints::getBarrierParam() const {
  return barrier_;
}


double Constraints::getFractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}

} // namespace robotoc