#ifndef ROBOTOC_IMPULSE_SPLIT_OCP_HXX_
#define ROBOTOC_IMPULSE_SPLIT_OCP_HXX_

#include "robotoc/impulse/impulse_split_ocp.hpp"

#include <cassert>

namespace robotoc {

inline ImpulseSplitOCP::ImpulseSplitOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, -1)),
    state_equation_(robot),
    impulse_dynamics_(robot),
    stage_cost_(0) {
}


inline ImpulseSplitOCP::ImpulseSplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    state_equation_(),
    impulse_dynamics_(),
    stage_cost_(0) {
}


inline ImpulseSplitOCP::~ImpulseSplitOCP() {
}


inline bool ImpulseSplitOCP::isFeasible(Robot& robot, 
                                        const ImpulseStatus& impulse_status,
                                        const ImpulseSplitSolution& s) {
  return constraints_->isFeasible(robot, impulse_status, constraints_data_, s);
}


inline void ImpulseSplitOCP::initConstraints(Robot& robot,
                                             const ImpulseStatus& impulse_status,
                                             const ImpulseSplitSolution& s) { 
  constraints_data_ = constraints_->createConstraintsData(robot, -1);
  constraints_->setSlackAndDual(robot, impulse_status, constraints_data_, s);
}


inline void ImpulseSplitOCP::initConstraints(const ImpulseSplitOCP& other) { 
  constraints_data_.copySlackAndDual(other.constraintsData());
}


inline const ConstraintsData& ImpulseSplitOCP::constraintsData() const {
  return constraints_data_;
}


inline void ImpulseSplitOCP::evalOCP(Robot& robot, 
                                     const ImpulseStatus& impulse_status, 
                                     const GridInfo& grid_info, 
                                     const ImpulseSplitSolution& s, 
                                     const Eigen::VectorXd& q_next, 
                                     const Eigen::VectorXd& v_next, 
                                     ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v+s.dv);
  stage_cost_ = cost_->evalImpulseCost(robot, impulse_status, cost_data_, 
                                       grid_info, s);
  constraints_->evalConstraint(robot, impulse_status, constraints_data_, s);
  stage_cost_ += constraints_data_.logBarrier();
  state_equation_.evalStateEquation(robot, s, q_next, v_next, kkt_residual);
  impulse_dynamics_.evalImpulseDynamics(robot, impulse_status, s);
}


inline void ImpulseSplitOCP::computeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const GridInfo& grid_info, 
    const Eigen::VectorXd& q_prev, const ImpulseSplitSolution& s, 
    const SplitSolution& s_next, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeImpulseCost(robot, impulse_status, cost_data_, 
                                            grid_info, s, kkt_residual);
  constraints_->linearizeConstraints(robot, impulse_status, constraints_data_, 
                                     s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  state_equation_.linearizeStateEquation(robot, q_prev, s, s_next, 
                                         kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s, 
                                             kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
}


inline void ImpulseSplitOCP::computeKKTSystem(
    Robot& robot, const ImpulseStatus& impulse_status, const GridInfo& grid_info,
    const Eigen::VectorXd& q_prev, const ImpulseSplitSolution& s, 
    const SplitSolution& s_next, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeImpulseCost(robot, impulse_status, cost_data_,  
                                             grid_info, s, kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, impulse_status, constraints_data_, 
                                     s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  state_equation_.linearizeStateEquation(robot, q_prev, s, s_next, 
                                         kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s, 
                                             kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
  constraints_->condenseSlackAndDual(impulse_status, constraints_data_, 
                                     kkt_matrix, kkt_residual);
  impulse_dynamics_.condenseImpulseDynamics(robot, impulse_status,
                                            kkt_matrix, kkt_residual);
  state_equation_.correctLinearizedStateEquation(robot, s, s_next, 
                                                 kkt_matrix, kkt_residual);
}


inline void ImpulseSplitOCP::expandPrimal(const ImpulseStatus& impulse_status, 
                                          ImpulseSplitDirection& d) {
  d.setImpulseStatus(impulse_status);
  impulse_dynamics_.expandPrimal(d);
  constraints_->expandSlackAndDual(impulse_status, constraints_data_, d);
}


inline void ImpulseSplitOCP::expandDual(const SplitDirection& d_next, 
                                        ImpulseSplitDirection& d) {
  impulse_dynamics_.expandDual(d_next, d);
  state_equation_.correctCostateDirection(d);
}


inline double ImpulseSplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double ImpulseSplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void ImpulseSplitOCP::updatePrimal(
    const Robot& robot, const double primal_step_size, 
    const ImpulseSplitDirection& d, ImpulseSplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void ImpulseSplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline double ImpulseSplitOCP::KKTError(
    const ImpulseSplitKKTResidual& kkt_residual) const {
  double err = 0;
  err += kkt_residual.KKTError();
  err += impulse_dynamics_.KKTError();
  err += constraints_data_.KKTError();
  return err;
}


inline double ImpulseSplitOCP::stageCost() const {
  return stage_cost_;
}


inline double ImpulseSplitOCP::constraintViolation(
    const ImpulseSplitKKTResidual& kkt_residual) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += constraints_data_.constraintViolation();
  vio += impulse_dynamics_.constraintViolation();
  return vio;
}

} // namespace robotoc

#endif // ROBOTOC_IMPULSE_SPLIT_OCP_HXX_ 