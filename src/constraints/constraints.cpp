#include "robotoc/constraints/constraints.hpp"
#include "robotoc/constraints/constraints_impl.hpp"

#include <stdexcept>
#include <cassert>


namespace robotoc {

Constraints::Constraints(const double barrier, 
                         const double fraction_to_boundary_rule) 
  : position_level_constraints_(),
    velocity_level_constraints_(),
    acceleration_level_constraints_(),
    impulse_level_constraints_(),
    barrier_(barrier), 
    fraction_to_boundary_rule_(fraction_to_boundary_rule) {
  if (barrier <= 0) {
    throw std::out_of_range(
        "[Constraints] invalid argment: 'barrier' must be positive!");
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


Constraints::~Constraints() {
}


void Constraints::push_back(ConstraintComponentBasePtr constraint_component) {
  constraint_component->setBarrier(barrier_);
  constraint_component->setFractionToBoundaryRule(fraction_to_boundary_rule_);
  if (constraint_component->kinematicsLevel() 
        == KinematicsLevel::PositionLevel) {
    position_level_constraints_.push_back(constraint_component);
  }
  else if (constraint_component->kinematicsLevel() 
              == KinematicsLevel::VelocityLevel) {
    velocity_level_constraints_.push_back(constraint_component);
  }
  else if (constraint_component->kinematicsLevel() 
              == KinematicsLevel::AccelerationLevel) {
    acceleration_level_constraints_.push_back(constraint_component);
  }
}


void Constraints::push_back(ImpulseConstraintComponentBasePtr constraint_component) {
  constraint_component->setBarrier(barrier_);
  constraint_component->setFractionToBoundaryRule(fraction_to_boundary_rule_);
  if (constraint_component->kinematicsLevel() 
        == KinematicsLevel::AccelerationLevel) {
    // Only the acceleration level constraints are valid at impulse stage.
    impulse_level_constraints_.push_back(constraint_component); 
  }
}


void Constraints::clear() {
  constraintsimpl::clear(position_level_constraints_);
  constraintsimpl::clear(velocity_level_constraints_);
  constraintsimpl::clear(acceleration_level_constraints_);
  constraintsimpl::clear(impulse_level_constraints_);
}


bool Constraints::useKinematics() const {
  if (constraintsimpl::useKinematics(position_level_constraints_)) {
    return true;
  }
  if (constraintsimpl::useKinematics(velocity_level_constraints_)) {
    return true;
  }
  if (constraintsimpl::useKinematics(acceleration_level_constraints_)) {
    return true;
  }
  return false;
}


ConstraintsData Constraints::createConstraintsData(const Robot& robot, 
                                                   const int time_stage) const {
  ConstraintsData data(time_stage);
  if (data.isPositionLevelValid()) {
    constraintsimpl::createConstraintsData(position_level_constraints_, 
                                           data.position_level_data);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::createConstraintsData(velocity_level_constraints_, 
                                           data.velocity_level_data);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::createConstraintsData(acceleration_level_constraints_, 
                                           data.acceleration_level_data);
  }
  if (data.isImpulseLevelValid()) {
    constraintsimpl::createConstraintsData(impulse_level_constraints_, 
                                           data.impulse_level_data);
  }
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


bool Constraints::isFeasible(Robot& robot, const ImpulseStatus& impulse_status, 
                             ConstraintsData& data, 
                             const ImpulseSplitSolution& s) const {
  if (data.isImpulseLevelValid()) {
    if (!constraintsimpl::isFeasible(impulse_level_constraints_, robot, 
                                     impulse_status, 
                                     data.impulse_level_data, s)) {
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
                                  const ImpulseStatus& impulse_status,
                                  ConstraintsData& data, 
                                  const ImpulseSplitSolution& s) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::setSlackAndDual(impulse_level_constraints_, robot,
                                     impulse_status, 
                                     data.impulse_level_data, s);
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
                                 const ImpulseStatus& impulse_status, 
                                 ConstraintsData& data, 
                                 const ImpulseSplitSolution& s) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::evalConstraint(impulse_level_constraints_, robot, 
                                    impulse_status, data.impulse_level_data, s);
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
                                       const ImpulseStatus& impulse_status, 
                                       ConstraintsData& data, 
                                       const ImpulseSplitSolution& s, 
                                       ImpulseSplitKKTResidual& kkt_residual) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::linearizeConstraints(impulse_level_constraints_, robot, 
                                          impulse_status, 
                                          data.impulse_level_data, s, kkt_residual);
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


void Constraints::condenseSlackAndDual(const ImpulseStatus& impulse_status, 
                                       ConstraintsData& data, 
                                       ImpulseSplitKKTMatrix& kkt_matrix, 
                                       ImpulseSplitKKTResidual& kkt_residual) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::condenseSlackAndDual(impulse_level_constraints_, 
                                          impulse_status, 
                                          data.impulse_level_data, 
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


void Constraints::expandSlackAndDual(const ImpulseStatus& impulse_status, 
                                     ConstraintsData& data, 
                                     const ImpulseSplitDirection& d) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::expandSlackAndDual(impulse_level_constraints_, 
                                        impulse_status, 
                                        data.impulse_level_data, d);
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
  if (data.isImpulseLevelValid()) {
    const double step_size = constraintsimpl::maxSlackStepSize(
        impulse_level_constraints_, data.impulse_level_data);
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
  if (data.isImpulseLevelValid()) {
    const double step_size = constraintsimpl::maxDualStepSize(
        impulse_level_constraints_, data.impulse_level_data);
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
  if (data.isImpulseLevelValid()) {
    constraintsimpl::updateSlack(data.impulse_level_data, step_size);
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
  if (data.isImpulseLevelValid()) {
    constraintsimpl::updateDual(data.impulse_level_data, step_size);
  }
}


void Constraints::setBarrier(const double _barrier) {
  assert(_barrier > 0);
  constraintsimpl::setBarrier(position_level_constraints_, _barrier);
  constraintsimpl::setBarrier(velocity_level_constraints_, _barrier);
  constraintsimpl::setBarrier(acceleration_level_constraints_, _barrier);
  constraintsimpl::setBarrier(impulse_level_constraints_, _barrier);
  barrier_ = _barrier;
}


void Constraints::setFractionToBoundaryRule(
    const double fraction_to_boundary_rule) {
  assert(fraction_to_boundary_rule > 0);
  constraintsimpl::setFractionToBoundaryRule(position_level_constraints_, 
                                             fraction_to_boundary_rule);
  constraintsimpl::setFractionToBoundaryRule(velocity_level_constraints_, 
                                             fraction_to_boundary_rule);
  constraintsimpl::setFractionToBoundaryRule(acceleration_level_constraints_, 
                                             fraction_to_boundary_rule);
  constraintsimpl::setFractionToBoundaryRule(impulse_level_constraints_, 
                                             fraction_to_boundary_rule);
  fraction_to_boundary_rule_ = fraction_to_boundary_rule;
}


double Constraints::barrier() const {
  return barrier_;
}


double Constraints::fractionToBoundaryRule() const {
  return fraction_to_boundary_rule_;
}

} // namespace robotoc