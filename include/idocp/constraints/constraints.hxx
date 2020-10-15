#ifndef IDOCP_CONSTRAINTS_HXX_
#define IDOCP_CONSTRAINTS_HXX_

#include <assert.h>

namespace idocp {

inline Constraints::Constraints() 
  : position_level_constraints_(),
    velocity_level_constraints_(),
    acceleration_level_constraints_() {
}


inline Constraints::~Constraints() {
}


inline void Constraints::push_back(
    const std::shared_ptr<ConstraintComponentBase>& constraint) {
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


inline void Constraints::clear() {
  clear_impl(position_level_constraints_);
  clear_impl(velocity_level_constraints_);
  clear_impl(acceleration_level_constraints_);
}


inline void Constraints::clear_impl(
    std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints) {
  constraints.clear();
}


inline bool Constraints::isEmpty() const {
  if (isEmpty_impl(position_level_constraints_)) {
    return true;
  }
  if (isEmpty_impl(velocity_level_constraints_)) {
    return true;
  }
  if (isEmpty_impl(acceleration_level_constraints_)) {
    return true;
  }
  return false;
}


inline bool Constraints::isEmpty_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints) {
  if (constraints.empty()) {
    return true;
  }
  else {
    return false;
  }
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints) {
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolution& s) {
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
                                         const double dtau, 
                                         const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    setSlackAndDual_impl(position_level_constraints_, robot, 
                         data.position_level_data, dtau, s);
  }
  if (data.isVelocityLevelValid()) {
    setSlackAndDual_impl(velocity_level_constraints_, robot, 
                         data.velocity_level_data, dtau, s);
  }
  setSlackAndDual_impl(acceleration_level_constraints_, robot,
                       data.acceleration_level_data, dtau, s);
}


inline void Constraints::setSlackAndDual_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->setSlackAndDual(robot, data[i], dtau, s);
  }
}


inline void Constraints::augmentDualResidual(Robot& robot, 
                                             ConstraintsData& data, 
                                             const double dtau, 
                                             const SplitSolution& s,
                                             KKTResidual& kkt_residual) const {
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s, KKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->augmentDualResidual(robot, data[i], dtau, s, kkt_residual);
  }
}


inline void Constraints::augmentDualResidual(const Robot& robot, 
                                             ConstraintsData& data, 
                                             const double dtau, 
                                             const Eigen::VectorXd& u,
                                             Eigen::VectorXd& lu) const {
  assert(data.acceleration_level_data.size() 
          == acceleration_level_constraints_.size());
  assert(u.size() == robot.dimv());
  assert(lu.size() == robot.dimv());
  for (int i=0; i<acceleration_level_constraints_.size(); ++i) {
    assert(data.acceleration_level_data[i].dimc() 
            == acceleration_level_constraints_[i]->dimc());
    assert(data.acceleration_level_data[i].checkDimensionalConsistency());
    acceleration_level_constraints_[i]->augmentDualResidual(
          robot, data.acceleration_level_data[i], dtau, u, lu);
  }
}


inline void Constraints::condenseSlackAndDual(Robot& robot,
                                              ConstraintsData& data, 
                                              const double dtau, 
                                              const SplitSolution& s,
                                              KKTMatrix& kkt_matrix, 
                                              KKTResidual& kkt_residual) const {
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->condenseSlackAndDual(robot, data[i], dtau, s, kkt_matrix, 
                                         kkt_residual);
  }
}


inline void Constraints::condenseSlackAndDual(const Robot& robot,
                                              ConstraintsData& data, 
                                              const double dtau, 
                                              const Eigen::VectorXd& u,
                                              Eigen::MatrixXd& Quu, 
                                              Eigen::VectorXd& lu) const {
  assert(data.acceleration_level_data.size() 
          == acceleration_level_constraints_.size());
  assert(u.size() == robot.dimv());
  assert(Quu.rows() == robot.dimv());
  assert(Quu.cols() == robot.dimv());
  assert(lu.size() == robot.dimv());
  for (int i=0; i<acceleration_level_constraints_.size(); ++i) {
    assert(data.acceleration_level_data[i].dimc() 
            == acceleration_level_constraints_[i]->dimc());
    assert(data.acceleration_level_data[i].checkDimensionalConsistency());
    acceleration_level_constraints_[i]->condenseSlackAndDual(
        robot, data.acceleration_level_data[i], dtau, u, Quu, lu);
  }
}


inline void Constraints::computeSlackAndDualDirection(
    Robot& robot, ConstraintsData& data, const double dtau, 
    const SplitSolution& s, const SplitDirection& d) const {
  if (data.isPositionLevelValid()) {
    computeSlackAndDualDirection_impl(position_level_constraints_, robot, 
                                      data.position_level_data, dtau, s, d);
  }
  if (data.isVelocityLevelValid()) {
    computeSlackAndDualDirection_impl(velocity_level_constraints_, robot, 
                                      data.velocity_level_data, dtau, s, d);
  }
  computeSlackAndDualDirection_impl(acceleration_level_constraints_, robot, 
                                    data.acceleration_level_data, dtau, s, d);
}


inline void Constraints::computeSlackAndDualDirection_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s, const SplitDirection& d) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computeSlackAndDualDirection(robot, data[i], dtau, s, d);
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    std::vector<ConstraintComponentData>& data, const double step_size) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->updateDual(data[i], step_size);
  }
}


inline double Constraints::costSlackBarrier(const ConstraintsData& data) const {
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
  return cost;
}


inline double Constraints::costSlackBarrier_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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
                                            const double step_size) const {
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
  return cost;
}


inline double Constraints::costSlackBarrier_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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
    Robot& robot, ConstraintsData& data, const double dtau, 
    const SplitSolution& s) const {
  if (data.isPositionLevelValid()) {
    computePrimalAndDualResidual_impl(position_level_constraints_, robot, 
                                      data.position_level_data, dtau, s);
  }
  if (data.isVelocityLevelValid()) {
    computePrimalAndDualResidual_impl(velocity_level_constraints_, robot, 
                                      data.velocity_level_data, dtau, s);
  }
  computePrimalAndDualResidual_impl(acceleration_level_constraints_, robot, 
                                    data.acceleration_level_data, dtau, s);
}


inline void Constraints::computePrimalAndDualResidual_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computePrimalAndDualResidual(robot, data[i], dtau, s);
  }
}


inline double Constraints::l1NormPrimalResidual(
    const ConstraintsData& data) const {
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
  return l1_norm;
}


inline double Constraints::l1NormPrimalResidual_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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
    const ConstraintsData& data) const {
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
  return squared_norm;
}


inline double Constraints::squaredNormPrimalAndDualResidual_impl(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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