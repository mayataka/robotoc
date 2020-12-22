#ifndef IDOCP_CONSTRAINTS_HXX_
#define IDOCP_CONSTRAINTS_HXX_

#include <cassert>

namespace idocp {

inline Constraints::Constraints() 
  : position_level_constraints_(),
    velocity_level_constraints_(),
    acceleration_level_constraints_(),
    impulse_constraints_(std::make_shared<ImpulseConstraints>()) {
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
  impulse_constraints_->push_back(constraint);
}


inline std::shared_ptr<ImpulseConstraints> 
Constraints::getImpulseConstraints() {
  return impulse_constraints_;
}


inline void Constraints::clear() {
  clear_impl(position_level_constraints_);
  clear_impl(velocity_level_constraints_);
  clear_impl(acceleration_level_constraints_);
  impulse_constraints_->clear();
}


inline void Constraints::clear_impl(
    std::vector<ConstraintComponentBasePtr>& constraints) {
  constraints.clear();
}


inline bool Constraints::useKinematics() const {
  if (useKinematics_impl(position_level_constraints_)) {
    return true;
  }
  if (useKinematics_impl(velocity_level_constraints_)) {
    return true;
  }
  if (useKinematics_impl(acceleration_level_constraints_)) {
    return true;
  }
  return false;
}


inline bool Constraints::useKinematics_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints) {
  for (const auto constraint : constraints) {
    if (constraint->useKinematics()) {
      return true;
    }
  }
  return false;
}


inline ConstraintsData Constraints::createConstraintsData(
    const Robot& robot, const int time_stage) const {
  assert(time_stage >= 0);
  ConstraintsData data(time_stage);
  if (data.isPositionLevelValid()) {
    for (const auto constraint : position_level_constraints_) {
      data.position_level_data.push_back(ConstraintComponentData(
          constraint->dimc()));
    }
  }
  if (data.isVelocityLevelValid()) {
    for (const auto constraint : velocity_level_constraints_) {
      data.velocity_level_data.push_back(ConstraintComponentData(
          constraint->dimc()));
    }
  }
  for (const auto constraint : acceleration_level_constraints_) {
    data.acceleration_level_data.push_back(ConstraintComponentData(
        constraint->dimc()));
  }
  return data;
}


inline bool Constraints::isFeasible(Robot& robot, ConstraintsData& data, 
                                    const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    if (!isFeasible_impl(position_level_constraints_, robot, 
                         data.position_level_data, s)) {
      return false;
    }
  }
  if (data.isVelocityLevelValid()) {
    if (!isFeasible_impl(velocity_level_constraints_, robot, 
                         data.velocity_level_data, s)) {
      return false;
    }
  }
  if (!isFeasible_impl(acceleration_level_constraints_, robot, 
                       data.acceleration_level_data, s)) {
    return false;
  }
  return true;
}


inline bool Constraints::isFeasible_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const SplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    bool feasible = constraints[i]->isFeasible(robot, data[i], s);
    if (!feasible) {
      return false;
    }
  }
  return true;
}


inline void Constraints::setSlackAndDual(Robot& robot, ConstraintsData& data, 
                                         const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    setSlackAndDual_impl(position_level_constraints_, robot, 
                         data.position_level_data, s);
  }
  if (data.isVelocityLevelValid()) {
    setSlackAndDual_impl(velocity_level_constraints_, robot, 
                         data.velocity_level_data, s);
  }
  setSlackAndDual_impl(acceleration_level_constraints_, robot,
                       data.acceleration_level_data, s);
}


inline void Constraints::setSlackAndDual_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const SplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->setSlackAndDual(robot, data[i], s);
  }
}


inline void Constraints::augmentDualResidual(
    Robot& robot, ConstraintsData& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (data.isPositionLevelValid()) {
    augmentDualResidual_impl(position_level_constraints_, robot, 
                             data.position_level_data, dtau, s, kkt_residual);
  }
  if (data.isVelocityLevelValid()) {
    augmentDualResidual_impl(velocity_level_constraints_, robot, 
                             data.velocity_level_data, dtau, s, kkt_residual);
  }
  augmentDualResidual_impl(acceleration_level_constraints_, robot, 
                           data.acceleration_level_data, dtau, s, kkt_residual);
}


inline void Constraints::augmentDualResidual_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->augmentDualResidual(robot, data[i], dtau, s, kkt_residual);
  }
}


inline void Constraints::condenseSlackAndDual(Robot& robot,
                                              ConstraintsData& data, 
                                              const double dtau, 
                                              const SplitSolution& s,
                                              SplitKKTMatrix& kkt_matrix, 
                                              SplitKKTResidual& kkt_residual) const {
  if (data.isPositionLevelValid()) {
    condenseSlackAndDual_impl(position_level_constraints_, robot,
                              data.position_level_data, dtau, s, kkt_matrix,
                              kkt_residual);
  }
  if (data.isVelocityLevelValid()) {
    condenseSlackAndDual_impl(velocity_level_constraints_, robot,
                              data.velocity_level_data, dtau, s, kkt_matrix,
                              kkt_residual);
  }
  condenseSlackAndDual_impl(acceleration_level_constraints_, robot,
                            data.acceleration_level_data, dtau, s, kkt_matrix,
                            kkt_residual);
}


inline void Constraints::condenseSlackAndDual_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->condenseSlackAndDual(robot, data[i], dtau, s, kkt_matrix, 
                                         kkt_residual);
  }
}


inline void Constraints::computeSlackAndDualDirection(
    Robot& robot, ConstraintsData& data, const SplitSolution& s, 
    const SplitDirection& d) const {
  if (data.isPositionLevelValid()) {
    computeSlackAndDualDirection_impl(position_level_constraints_, robot, 
                                      data.position_level_data, s, d);
  }
  if (data.isVelocityLevelValid()) {
    computeSlackAndDualDirection_impl(velocity_level_constraints_, robot, 
                                      data.velocity_level_data, s, d);
  }
  computeSlackAndDualDirection_impl(acceleration_level_constraints_, robot, 
                                    data.acceleration_level_data, s, d);
}


inline void Constraints::computeSlackAndDualDirection_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const SplitSolution& s, 
    const SplitDirection& d) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computeSlackAndDualDirection(robot, data[i], s, d);
  }
}


inline double Constraints::maxSlackStepSize(const ConstraintsData& data) const {
  double min_step_size = 1;
  if (data.isPositionLevelValid()) {
    const double step_size = maxSlackStepSize_impl(position_level_constraints_,
                                                   data.position_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isVelocityLevelValid()) {
    const double step_size = maxSlackStepSize_impl(velocity_level_constraints_,
                                                   data.velocity_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  const double step_size = maxSlackStepSize_impl(acceleration_level_constraints_,
                                                 data.acceleration_level_data);
  if (step_size < min_step_size) {
    min_step_size = step_size;
  }
  return min_step_size;
}


inline double Constraints::maxSlackStepSize_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double min_step_size = 1;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    const double step_size = constraints[i]->maxSlackStepSize(data[i]);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


inline double Constraints::maxDualStepSize(const ConstraintsData& data) const {
  double min_step_size = 1;
  if (data.isPositionLevelValid()) {
    const double step_size = maxDualStepSize_impl(position_level_constraints_,
                                                  data.position_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  if (data.isVelocityLevelValid()) {
    const double step_size = maxDualStepSize_impl(velocity_level_constraints_,
                                                  data.velocity_level_data);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  const double step_size = maxDualStepSize_impl(acceleration_level_constraints_,
                                                data.acceleration_level_data);
  if (step_size < min_step_size) {
    min_step_size = step_size;
  }
  return min_step_size;
}


inline double Constraints::maxDualStepSize_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double min_step_size = 1;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    const double step_size = constraints[i]->maxDualStepSize(data[i]);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


inline void Constraints::updateSlack(ConstraintsData& data, 
                                     const double step_size) const {
  assert(step_size >= 0);
  assert(step_size <= 1);
  if (data.isPositionLevelValid()) {
    updateSlack_impl(position_level_constraints_, data.position_level_data, 
                     step_size);
  }
  if (data.isVelocityLevelValid()) {
    updateSlack_impl(velocity_level_constraints_, data.velocity_level_data, 
                     step_size);
  }
  updateSlack_impl(acceleration_level_constraints_, data.acceleration_level_data, 
                   step_size);
}


inline void Constraints::updateSlack_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    std::vector<ConstraintComponentData>& data, const double step_size) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->updateSlack(data[i], step_size);
  }
}


inline void Constraints::updateDual(ConstraintsData& data, 
                                    const double step_size) const {
  assert(step_size >= 0);
  assert(step_size <= 1);
  if (data.isPositionLevelValid()) {
    updateDual_impl(position_level_constraints_, data.position_level_data, 
                    step_size);
  }
  if (data.isVelocityLevelValid()) {
    updateDual_impl(velocity_level_constraints_, data.velocity_level_data, 
                    step_size);
  }
  updateDual_impl(acceleration_level_constraints_, data.acceleration_level_data, 
                  step_size);
}


inline void Constraints::updateDual_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    std::vector<ConstraintComponentData>& data, const double step_size) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->updateDual(data[i], step_size);
  }
}


inline double Constraints::costSlackBarrier(const ConstraintsData& data,
                                            const double dtau) const {
  assert(dtau >= 0);
  double cost = 0;
  if (data.isPositionLevelValid()) {
    cost += costSlackBarrier_impl(position_level_constraints_, 
                                  data.position_level_data);
  }
  if (data.isVelocityLevelValid()) {
    cost += costSlackBarrier_impl(velocity_level_constraints_, 
                                  data.velocity_level_data);
  }
  cost += costSlackBarrier_impl(acceleration_level_constraints_, 
                                data.acceleration_level_data);
  return dtau * cost;
}


inline double Constraints::costSlackBarrier_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double cost = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    cost += constraints[i]->costSlackBarrier(data[i]);
  }
  return cost;
}


inline double Constraints::costSlackBarrier(const ConstraintsData& data, 
                                            const double dtau,
                                            const double step_size) const {
  assert(dtau >= 0);
  assert(step_size >= 0);
  assert(step_size <= 1);
  double cost = 0;
  if (data.isPositionLevelValid()) {
    cost += costSlackBarrier_impl(position_level_constraints_, 
                                  data.position_level_data, step_size);
  }
  if (data.isVelocityLevelValid()) {
    cost += costSlackBarrier_impl(velocity_level_constraints_, 
                                  data.velocity_level_data, step_size);
  }
  cost += costSlackBarrier_impl(acceleration_level_constraints_, 
                                data.acceleration_level_data, step_size);
  return dtau * cost;
}


inline double Constraints::costSlackBarrier_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    const std::vector<ConstraintComponentData>& data, const double step_size) {
  assert(constraints.size() == data.size());
  double cost = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    cost += constraints[i]->costSlackBarrier(data[i], step_size);
  }
  return cost;
}


inline void Constraints::computePrimalAndDualResidual(
    Robot& robot, ConstraintsData& data, const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    computePrimalAndDualResidual_impl(position_level_constraints_, robot, 
                                      data.position_level_data, s);
  }
  if (data.isVelocityLevelValid()) {
    computePrimalAndDualResidual_impl(velocity_level_constraints_, robot, 
                                      data.velocity_level_data, s);
  }
  computePrimalAndDualResidual_impl(acceleration_level_constraints_, robot, 
                                    data.acceleration_level_data, s);
}


inline void Constraints::computePrimalAndDualResidual_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const SplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computePrimalAndDualResidual(robot, data[i], s);
  }
}


inline double Constraints::l1NormPrimalResidual(const ConstraintsData& data, 
                                                const double dtau) const {
  double l1_norm = 0;
  if (data.isPositionLevelValid()) {
    l1_norm += l1NormPrimalResidual_impl(position_level_constraints_, 
                                         data.position_level_data);
  }
  if (data.isVelocityLevelValid()) {
    l1_norm += l1NormPrimalResidual_impl(velocity_level_constraints_, 
                                         data.velocity_level_data);
  }
  l1_norm += l1NormPrimalResidual_impl(acceleration_level_constraints_, 
                                       data.acceleration_level_data);
  return dtau * l1_norm;
}


inline double Constraints::l1NormPrimalResidual_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double l1_norm = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    l1_norm += constraints[i]->l1NormPrimalResidual(data[i]);
  }
  return l1_norm;
}


inline double Constraints::squaredNormPrimalAndDualResidual(
    const ConstraintsData& data, const double dtau) const {
  double squared_norm = 0;
  if (data.isPositionLevelValid()) {
    squared_norm += squaredNormPrimalAndDualResidual_impl(
        position_level_constraints_, data.position_level_data);
  }
  if (data.isVelocityLevelValid()) {
    squared_norm += squaredNormPrimalAndDualResidual_impl(
        velocity_level_constraints_, data.velocity_level_data);
  }
  squared_norm += squaredNormPrimalAndDualResidual_impl(
      acceleration_level_constraints_, data.acceleration_level_data);
  return dtau * dtau * squared_norm;
}


inline double Constraints::squaredNormPrimalAndDualResidual_impl(
    const std::vector<ConstraintComponentBasePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double squared_norm = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    squared_norm += constraints[i]->squaredNormPrimalAndDualResidual(data[i]);
  }
  return squared_norm;
}

} // namespace idocp

#endif // IDOCP_CONSTRAINTS_HXX_