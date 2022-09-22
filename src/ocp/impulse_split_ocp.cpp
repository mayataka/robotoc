#include "robotoc/ocp/impulse_split_ocp.hpp"
#include "robotoc/dynamics/impulse_state_equation.hpp"
#include "robotoc/dynamics/impulse_dynamics.hpp"

#include <cassert>

namespace robotoc {

ImpulseSplitOCP::ImpulseSplitOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, -1)),
    state_equation_data_(robot),
    contact_dynamics_data_(robot),
    stage_cost_(0),
    barrier_cost_(0) {
}


ImpulseSplitOCP::ImpulseSplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    state_equation_data_(),
    contact_dynamics_data_(),
    stage_cost_(0),
    barrier_cost_(0) {
}


ImpulseSplitOCP::~ImpulseSplitOCP() {
}


bool ImpulseSplitOCP::isFeasible(Robot& robot, 
                                 const ImpulseStatus& impulse_status,
                                 const SplitSolution& s) {
  return constraints_->isFeasible(robot, impulse_status, constraints_data_, s);
}


void ImpulseSplitOCP::initConstraints(Robot& robot,
                                      const ImpulseStatus& impulse_status,
                                      const SplitSolution& s) { 
  constraints_data_ = constraints_->createConstraintsData(robot, -1);
  constraints_->setSlackAndDual(robot, impulse_status, constraints_data_, s);
}


void ImpulseSplitOCP::initConstraints(const ImpulseSplitOCP& other) { 
  constraints_data_.copySlackAndDual(other.constraintsData());
}


const ConstraintsData& ImpulseSplitOCP::constraintsData() const {
  return constraints_data_;
}


void ImpulseSplitOCP::evalOCP(Robot& robot, const ImpulseStatus& impulse_status, 
                              const GridInfo& grid_info, 
                              const SplitSolution& s, 
                              const Eigen::VectorXd& q_next, 
                              const Eigen::VectorXd& v_next, 
                              SplitKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  kkt_residual.setContactStatus(impulse_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v+s.dv);
  stage_cost_ = cost_->evalImpulseCost(robot, impulse_status, cost_data_, 
                                       grid_info, s);
  constraints_->evalConstraint(robot, impulse_status, constraints_data_, s);
  barrier_cost_ = constraints_data_.logBarrier();
  ImpulseStateEquation::eval(robot, s, q_next, v_next, kkt_residual);
  ImpulseDynamics::eval(robot, impulse_status, contact_dynamics_data_, s);
}


void ImpulseSplitOCP::computeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const GridInfo& grid_info, 
    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_matrix.setContactStatus(impulse_status);
  kkt_residual.setContactStatus(impulse_status);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeImpulseCost(robot, impulse_status, cost_data_, 
                                            grid_info, s, kkt_residual);
  constraints_->linearizeConstraints(robot, impulse_status, constraints_data_, 
                                     s, kkt_residual);
  barrier_cost_ = constraints_data_.logBarrier();
  ImpulseStateEquation::linearize(robot, state_equation_data_, q_prev, 
                                  s, s_next, kkt_matrix, kkt_residual);
  ImpulseDynamics::linearize(robot, impulse_status, contact_dynamics_data_, s, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
}


void ImpulseSplitOCP::computeKKTSystem(
    Robot& robot, const ImpulseStatus& impulse_status, const GridInfo& grid_info,
    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_matrix.setContactStatus(impulse_status);
  kkt_residual.setContactStatus(impulse_status);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeImpulseCost(robot, impulse_status, cost_data_,  
                                             grid_info, s, kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, impulse_status, constraints_data_, 
                                     s, kkt_residual);
  barrier_cost_ = constraints_data_.logBarrier();
  ImpulseStateEquation::linearize(robot, state_equation_data_, q_prev, 
                                  s, s_next, kkt_matrix, kkt_residual);
  ImpulseDynamics::linearize(robot, impulse_status, contact_dynamics_data_, s, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
  constraints_->condenseSlackAndDual(impulse_status, constraints_data_, 
                                     kkt_matrix, kkt_residual);
  ImpulseDynamics::condense(robot, impulse_status, contact_dynamics_data_, kkt_matrix, kkt_residual);
  ImpulseStateEquation::correctLinearize(robot, state_equation_data_,  
                                         s, s_next, kkt_matrix, kkt_residual);
}


void ImpulseSplitOCP::expandPrimal(const ImpulseStatus& impulse_status, 
                                   SplitDirection& d) {
  d.setContactStatus(impulse_status);
  ImpulseDynamics::expandPrimal(contact_dynamics_data_, d);
  constraints_->expandSlackAndDual(impulse_status, constraints_data_, d);
}


void ImpulseSplitOCP::expandDual(const SplitDirection& d_next, 
                                 SplitDirection& d) {
  ImpulseDynamics::expandDual(contact_dynamics_data_, d_next, d);
  ImpulseStateEquation::correctCostateDirection(state_equation_data_, d);
}


double ImpulseSplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double ImpulseSplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


void ImpulseSplitOCP::updatePrimal(const Robot& robot, 
                                   const double primal_step_size, 
                                   const SplitDirection& d, 
                                   SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d, true);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


void ImpulseSplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


double ImpulseSplitOCP::KKTError(
    const SplitKKTResidual& kkt_residual) const {
  double err = 0;
  err += kkt_residual.KKTError();
  err += contact_dynamics_data_.KKTError();
  err += constraints_data_.KKTError();
  return err;
}


double ImpulseSplitOCP::stageCost(const bool include_cost_barrier) const {
  if (include_cost_barrier) {
    return stage_cost_ + barrier_cost_; 
  }
  else {
    return stage_cost_;
  }
}


double ImpulseSplitOCP::constraintViolation(
    const SplitKKTResidual& kkt_residual) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += constraints_data_.constraintViolation();
  vio += contact_dynamics_data_.constraintViolation();
  return vio;
}

} // namespace robotoc