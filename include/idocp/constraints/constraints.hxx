#ifndef IDOCP_CONSTRAINTS_HXX_
#define IDOCP_CONSTRAINTS_HXX_

#include <cassert>

#include "idocp/constraints/constraints_impl.hpp"


namespace idocp {

inline Constraints::Constraints() 
  : position_level_constraints_(),
    velocity_level_constraints_(),
    acceleration_level_constraints_(),
    impulse_level_constraints_() {
}


inline Constraints::~Constraints() {
}


inline void Constraints::push_back(
    const ConstraintComponentBasePtr& constraint) {
  if (constraint->kinematicsLevel() == KinematicsLevel::PositionLevel) {
    position_level_constraints_.push_back(constraint);
  }
  else if (constraint->kinematicsLevel() == KinematicsLevel::VelocityLevel) {
    velocity_level_constraints_.push_back(constraint);
  }
  else if (constraint->kinematicsLevel() == KinematicsLevel::AccelerationLevel) {
    acceleration_level_constraints_.push_back(constraint);
  }
}


inline void Constraints::push_back(
    const ImpulseConstraintComponentBasePtr& constraint) {
  if (constraint->kinematicsLevel() == KinematicsLevel::AccelerationLevel) {
    // Only the acceleration level constraints are valid at impulse stage.
    impulse_level_constraints_.push_back(constraint); 
  }
}


inline void Constraints::clear() {
  constraintsimpl::clear(position_level_constraints_);
  constraintsimpl::clear(velocity_level_constraints_);
  constraintsimpl::clear(acceleration_level_constraints_);
  constraintsimpl::clear(impulse_level_constraints_);
}


inline bool Constraints::useKinematics() const {
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


inline ConstraintsData Constraints::createConstraintsData(
    const Robot& robot, const int time_stage) const {
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


inline bool Constraints::isFeasible(Robot& robot, ConstraintsData& data, 
                                    const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    if (!constraintsimpl::isFeasible(position_level_constraints_, robot, 
                                     data.position_level_data, s)) {
      return false;
    }
  }
  if (data.isVelocityLevelValid()) {
    if (!constraintsimpl::isFeasible(velocity_level_constraints_, robot, 
                                     data.velocity_level_data, s)) {
      return false;
    }
  }
  if (data.isAccelerationLevelValid()) {
    if (!constraintsimpl::isFeasible(acceleration_level_constraints_, robot, 
                                     data.acceleration_level_data, s)) {
      return false;
    }
  }
  return true;
}


inline bool Constraints::isFeasible(Robot& robot, ConstraintsData& data, 
                                    const ImpulseSplitSolution& s) const {
  if (data.isImpulseLevelValid()) {
    if (!constraintsimpl::isFeasible(impulse_level_constraints_, robot, 
                                     data.impulse_level_data, s)) {
      return false;
    }
  }
  return true;
}


inline void Constraints::setSlackAndDual(Robot& robot, ConstraintsData& data, 
                                         const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::setSlackAndDual(position_level_constraints_, robot, 
                                     data.position_level_data, s);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::setSlackAndDual(velocity_level_constraints_, robot, 
                                     data.velocity_level_data, s);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::setSlackAndDual(acceleration_level_constraints_, robot,
                                     data.acceleration_level_data, s);
  }
}


inline void Constraints::setSlackAndDual(Robot& robot, ConstraintsData& data, 
                                         const ImpulseSplitSolution& s) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::setSlackAndDual(impulse_level_constraints_, robot,
                                     data.impulse_level_data, s);
  }
}


inline void Constraints::augmentDualResidual(
    Robot& robot, ConstraintsData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::augmentDualResidual(position_level_constraints_, robot, 
                                         data.position_level_data, dtau, s, 
                                         kkt_residual);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::augmentDualResidual(velocity_level_constraints_, robot, 
                                         data.velocity_level_data, dtau, s, 
                                         kkt_residual);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::augmentDualResidual(acceleration_level_constraints_, robot, 
                                         data.acceleration_level_data, dtau, s, 
                                         kkt_residual);
  }
}


inline void Constraints::augmentDualResidual(
    Robot& robot, ConstraintsData& data, const ImpulseSplitSolution& s,
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::augmentDualResidual(impulse_level_constraints_, robot, 
                                         data.impulse_level_data, s, 
                                         kkt_residual);
  }
}


inline void Constraints::condenseSlackAndDual(
    Robot& robot, ConstraintsData& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::condenseSlackAndDual(position_level_constraints_, robot,
                                          data.position_level_data, dtau, s, 
                                          kkt_matrix, kkt_residual);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::condenseSlackAndDual(velocity_level_constraints_, robot,
                                          data.velocity_level_data, dtau, s, 
                                          kkt_matrix, kkt_residual);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::condenseSlackAndDual(acceleration_level_constraints_, robot,
                                          data.acceleration_level_data, dtau, s, 
                                          kkt_matrix, kkt_residual);
  }
}


inline void Constraints::condenseSlackAndDual(
    Robot& robot, ConstraintsData& data, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::condenseSlackAndDual(impulse_level_constraints_, robot,
                                          data.impulse_level_data, s, 
                                          kkt_matrix, kkt_residual);
  }
}


inline void Constraints::computeSlackAndDualDirection(
    Robot& robot, ConstraintsData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::computeSlackAndDualDirection(
        position_level_constraints_, robot, data.position_level_data, s, d);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::computeSlackAndDualDirection(
        velocity_level_constraints_, robot, data.velocity_level_data, s, d);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::computeSlackAndDualDirection(
        acceleration_level_constraints_, robot, data.acceleration_level_data, 
        s, d);
  }
}


inline void Constraints::computeSlackAndDualDirection(
    Robot& robot, ConstraintsData& data, const ImpulseSplitSolution& s, 
    const ImpulseSplitDirection& d) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::computeSlackAndDualDirection(
        impulse_level_constraints_, robot, data.impulse_level_data, s, d);
  }
}


inline double Constraints::maxSlackStepSize(const ConstraintsData& data) const {
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


inline double Constraints::maxDualStepSize(const ConstraintsData& data) const {
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


inline void Constraints::updateSlack(ConstraintsData& data, 
                                     const double step_size) {
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


inline void Constraints::updateDual(ConstraintsData& data, 
                                    const double step_size) {
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


inline double Constraints::costSlackBarrier(const ConstraintsData& data) const {
  double cost = 0;
  if (data.isPositionLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(position_level_constraints_, 
                                              data.position_level_data);
  }
  if (data.isVelocityLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(velocity_level_constraints_, 
                                              data.velocity_level_data);
  }
  if (data.isAccelerationLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(acceleration_level_constraints_, 
                                              data.acceleration_level_data);
  }
  if (data.isImpulseLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(impulse_level_constraints_, 
                                              data.impulse_level_data);
  }
  return cost;
}


inline double Constraints::costSlackBarrier(const ConstraintsData& data, 
                                            const double step_size) const {
  assert(step_size >= 0);
  assert(step_size <= 1);
  double cost = 0;
  if (data.isPositionLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(position_level_constraints_, 
                                              data.position_level_data, 
                                              step_size);
  }
  if (data.isVelocityLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(velocity_level_constraints_, 
                                              data.velocity_level_data, 
                                              step_size);
  }
  if (data.isAccelerationLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(acceleration_level_constraints_, 
                                              data.acceleration_level_data, 
                                              step_size);
  }
  if (data.isImpulseLevelValid()) {
    cost += constraintsimpl::costSlackBarrier(impulse_level_constraints_, 
                                              data.impulse_level_data, 
                                              step_size);
  }
  return cost;
}


inline void Constraints::computePrimalAndDualResidual(
    Robot& robot, ConstraintsData& data, const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    constraintsimpl::computePrimalAndDualResidual(
        position_level_constraints_, robot, data.position_level_data, s);
  }
  if (data.isVelocityLevelValid()) {
    constraintsimpl::computePrimalAndDualResidual(
        velocity_level_constraints_, robot, data.velocity_level_data, s);
  }
  if (data.isAccelerationLevelValid()) {
    constraintsimpl::computePrimalAndDualResidual(
        acceleration_level_constraints_, robot, data.acceleration_level_data, s);
  }
}


inline void Constraints::computePrimalAndDualResidual(
    Robot& robot, ConstraintsData& data, const ImpulseSplitSolution& s) const {
  if (data.isImpulseLevelValid()) {
    constraintsimpl::computePrimalAndDualResidual(
        impulse_level_constraints_, robot, data.impulse_level_data, s);
  }
}


inline double Constraints::l1NormPrimalResidual(
    const ConstraintsData& data) {
  double l1_norm = 0;
  if (data.isPositionLevelValid()) {
    l1_norm += constraintsimpl::l1NormPrimalResidual(data.position_level_data);
  }
  if (data.isVelocityLevelValid()) {
    l1_norm += constraintsimpl::l1NormPrimalResidual(data.velocity_level_data);
  }
  if (data.isAccelerationLevelValid()) {
    l1_norm += constraintsimpl::l1NormPrimalResidual(data.acceleration_level_data);
  }
  if (data.isImpulseLevelValid()) {
    l1_norm += constraintsimpl::l1NormPrimalResidual(data.impulse_level_data);
  }
  return l1_norm;
}


inline double Constraints::squaredNormPrimalAndDualResidual(
    const ConstraintsData& data) {
  double squared_norm = 0;
  if (data.isPositionLevelValid()) {
    squared_norm += constraintsimpl::squaredNormPrimalAndDualResidual(
        data.position_level_data);
  }
  if (data.isVelocityLevelValid()) {
    squared_norm += constraintsimpl::squaredNormPrimalAndDualResidual(
        data.velocity_level_data);
  }
  if (data.isAccelerationLevelValid()) {
    squared_norm += constraintsimpl::squaredNormPrimalAndDualResidual(
        data.acceleration_level_data);
  }
  if (data.isImpulseLevelValid()) {
    squared_norm += constraintsimpl::squaredNormPrimalAndDualResidual(
        data.impulse_level_data);
  }
  return squared_norm;
}


inline void Constraints::setBarrier(const double barrier) {
  constraintsimpl::setBarrier(position_level_constraints_, barrier);
  constraintsimpl::setBarrier(velocity_level_constraints_, barrier);
  constraintsimpl::setBarrier(acceleration_level_constraints_, barrier);
  constraintsimpl::setBarrier(impulse_level_constraints_, barrier);
}


inline void Constraints::setFractionToBoundaryRate(
    const double fraction_to_boundary_rate) {
  constraintsimpl::setFractionToBoundaryRate(position_level_constraints_, 
                                             fraction_to_boundary_rate);
  constraintsimpl::setFractionToBoundaryRate(velocity_level_constraints_, 
                                             fraction_to_boundary_rate);
  constraintsimpl::setFractionToBoundaryRate(acceleration_level_constraints_, 
                                             fraction_to_boundary_rate);
  constraintsimpl::setFractionToBoundaryRate(impulse_level_constraints_, 
                                             fraction_to_boundary_rate);
}

} // namespace idocp

#endif // IDOCP_CONSTRAINTS_HXX_