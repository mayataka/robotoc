#include "robotoc/unconstr/split_unconstr_ocp.hpp"

#include <stdexcept>
#include <cassert>


namespace robotoc {

SplitUnconstrOCP::SplitUnconstrOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    unconstr_dynamics_(robot),
    contact_status_(robot.createContactStatus()),
    stage_cost_(0) {
}


SplitUnconstrOCP::SplitUnconstrOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstr_dynamics_(),
    contact_status_(),
    stage_cost_(0) {
}


SplitUnconstrOCP::~SplitUnconstrOCP() {
}


bool SplitUnconstrOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, contact_status_, constraints_data_, s);
}


void SplitUnconstrOCP::initConstraints(Robot& robot, const int time_step, 
                                       const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, contact_status_, constraints_data_, s);
}


void SplitUnconstrOCP::evalOCP(Robot& robot, const GridInfo& grid_info,
                               const SplitSolution& s, 
                               const Eigen::VectorXd& q_next, 
                               const Eigen::VectorXd& v_next, 
                               SplitKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  robot.updateKinematics(s.q);
  kkt_residual.setZero();
  stage_cost_ = cost_->evalStageCost(robot, contact_status_, cost_data_, 
                                     grid_info, s);
  constraints_->evalConstraint(robot, contact_status_, constraints_data_, s);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::computeForwardEulerResidual(grid_info.dt, s, q_next, 
                                                       v_next, kkt_residual);
  unconstr_dynamics_.evalUnconstrDynamics(robot, s);
}


void SplitUnconstrOCP::computeKKTResidual(Robot& robot, 
                                          const GridInfo& grid_info,
                                          const SplitSolution& s,
                                          const SplitSolution& s_next,
                                          SplitKKTMatrix& kkt_matrix, 
                                          SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q);
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeStageCost(robot, contact_status_, cost_data_, 
                                          grid_info, s, kkt_residual);
  constraints_->linearizeConstraints(robot, contact_status_, constraints_data_, s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::linearizeForwardEuler(grid_info.dt, s, s_next, 
                                                 kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual);
  // kkt_residual.kkt_error = KKTError(kkt_residual, grid_info.dt);
}


void SplitUnconstrOCP::computeKKTSystem(Robot& robot, const GridInfo& grid_info,
                                        const SplitSolution& s, 
                                        const SplitSolution& s_next, 
                                        SplitKKTMatrix& kkt_matrix,
                                        SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeStageCost(robot, contact_status_, cost_data_, 
                                           grid_info, s, kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, contact_status_, constraints_data_, 
                                     s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  unconstr::stateequation::linearizeForwardEuler(grid_info.dt, s, s_next, 
                                                 kkt_matrix, kkt_residual);
  unconstr_dynamics_.linearizeUnconstrDynamics(robot, grid_info.dt, s, kkt_residual);
  // kkt_residual.kkt_error = KKTError(kkt_residual, grid_info.dt);
  constraints_->condenseSlackAndDual(contact_status_, constraints_data_, 
                                     kkt_matrix, kkt_residual);
  unconstr_dynamics_.condenseUnconstrDynamics(kkt_matrix, kkt_residual);
}


void SplitUnconstrOCP::expandPrimalAndDual(const double dt, 
                                           const SplitSolution& s, 
                                           const SplitKKTMatrix& kkt_matrix, 
                                           const SplitKKTResidual& kkt_residual, 
                                           SplitDirection& d) {
  assert(dt > 0);
  unconstr_dynamics_.expandPrimal(d);
  unconstr_dynamics_.expandDual(dt, kkt_matrix, kkt_residual, d);
  constraints_->expandSlackAndDual(contact_status_, constraints_data_, d);
}


double SplitUnconstrOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitUnconstrOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


void SplitUnconstrOCP::updatePrimal(const Robot& robot, 
                                    const double primal_step_size, 
                                    const SplitDirection& d, SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


void SplitUnconstrOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


double SplitUnconstrOCP::KKTError(const SplitKKTResidual& kkt_residual, 
                                  const double dt) const {
  assert(dt > 0);
  double err = 0;
  err += kkt_residual.KKTError();
  err += (dt*dt) * unconstr_dynamics_.KKTError();
  err += constraints_data_.KKTError();
  return err;
}


double SplitUnconstrOCP::stageCost() const {
  return stage_cost_;
}


double SplitUnconstrOCP::constraintViolation(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += dt * unconstr_dynamics_.constraintViolation();
  vio += constraints_data_.constraintViolation();
  return vio;
}

} // namespace robotoc