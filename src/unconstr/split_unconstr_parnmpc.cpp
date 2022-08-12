#include "robotoc/unconstr/split_unconstr_parnmpc.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>


namespace robotoc {

SplitUnconstrParNMPC::SplitUnconstrParNMPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    unconstr_dynamics_(robot),
    contact_status_(robot.createContactStatus()),
    use_kinematics_(false),
    stage_cost_(0) {
  if (cost_->useKinematics() || constraints_->useKinematics()) {
    use_kinematics_ = true;
  }
  if (robot.hasFloatingBase()) {
    throw std::logic_error(
        "[SplitUnconstrParNMPC] robot has floating base: robot should have no constraints!");
  }
  if (robot.maxNumContacts() > 0) {
    throw std::logic_error(
        "[SplitUnconstrParNMPC] robot can have contacts: robot should have no constraints!");
  }
}


SplitUnconstrParNMPC::SplitUnconstrParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstr_dynamics_(),
    contact_status_(),
    use_kinematics_(false),
    stage_cost_(0) {
}


SplitUnconstrParNMPC::~SplitUnconstrParNMPC() {
}


bool SplitUnconstrParNMPC::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, contact_status_, constraints_data_, s);
}


void SplitUnconstrParNMPC::initConstraints(Robot& robot, const int time_step, 
                                           const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, contact_status_, constraints_data_, s);
}


void SplitUnconstrParNMPC::evalOCP(Robot& robot, const GridInfo& grid_info,
                                   const Eigen::VectorXd& q_prev, 
                                   const Eigen::VectorXd& v_prev, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_residual.setZero();
  stage_cost_ = cost_->evalStageCost(robot, contact_status_, cost_data_, 
                                     grid_info, s);
  constraints_->evalConstraint(robot, contact_status_, constraints_data_, s);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::computeBackwardEulerResidual(grid_info.dt, q_prev, v_prev, 
                                                        s, kkt_residual);
  unconstr_dynamics_.evalUnconstrDynamics(robot, s);
}


void SplitUnconstrParNMPC::computeKKTResidual(Robot& robot, 
                                              const GridInfo& grid_info,
                                              const Eigen::VectorXd& q_prev, 
                                              const Eigen::VectorXd& v_prev, 
                                              const SplitSolution& s, 
                                              const SplitSolution& s_next, 
                                              SplitKKTMatrix& kkt_matrix, 
                                              SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeStageCost(robot, contact_status_, cost_data_, 
                                          grid_info, s, kkt_residual);
  constraints_->linearizeConstraints(robot, contact_status_, constraints_data_, 
                                     s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::linearizeBackwardEuler(grid_info.dt, q_prev, v_prev, 
                                                  s, s_next, kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual, grid_info.dt);
}


void SplitUnconstrParNMPC::computeKKTSystem(Robot& robot,  
                                            const GridInfo& grid_info,
                                            const Eigen::VectorXd& q_prev, 
                                            const Eigen::VectorXd& v_prev, 
                                            const SplitSolution& s, 
                                            const SplitSolution& s_next, 
                                            SplitKKTMatrix& kkt_matrix,
                                            SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  if (use_kinematics_) {
    robot.updateKinematics(s.q);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeStageCost(robot, contact_status_, cost_data_, 
                                           grid_info, s, kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, contact_status_, constraints_data_, 
                                     s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::linearizeBackwardEuler(grid_info.dt, q_prev, v_prev, 
                                                  s, s_next, kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual, grid_info.dt);
  constraints_->condenseSlackAndDual(contact_status_, constraints_data_, 
                                     kkt_matrix, kkt_residual);
  unconstr_dynamics_.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


void SplitUnconstrParNMPC::expandPrimalAndDual(const double dt, 
                                               const SplitSolution& s, 
                                               const SplitKKTMatrix& kkt_matrix, 
                                               const SplitKKTResidual& kkt_residual, 
                                               SplitDirection& d) {
  assert(dt > 0);
  unconstr_dynamics_.expandPrimal(d);
  unconstr_dynamics_.expandDual(dt, kkt_matrix, kkt_residual, d);
  constraints_->expandSlackAndDual(contact_status_, constraints_data_, d);
}


double SplitUnconstrParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitUnconstrParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


void SplitUnconstrParNMPC::updatePrimal(const Robot& robot, 
                                        const double primal_step_size, 
                                        const SplitDirection& d, 
                                        SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


void SplitUnconstrParNMPC::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


double SplitUnconstrParNMPC::KKTError(const SplitKKTResidual& kkt_residual, 
                                      const double dt) const {
  assert(dt > 0);
  double err = 0;
  err += kkt_residual.KKTError();
  err += (dt*dt) * unconstr_dynamics_.KKTError();
  err += constraints_data_.KKTError();
  return err;
}


double SplitUnconstrParNMPC::stageCost() const {
  return stage_cost_;
}


double SplitUnconstrParNMPC::constraintViolation(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += dt * unconstr_dynamics_.constraintViolation();
  vio += constraints_data_.constraintViolation();
  return vio;
}

} // namespace robotoc