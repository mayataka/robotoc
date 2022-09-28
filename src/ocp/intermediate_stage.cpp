#include "robotoc/ocp/intermediate_stage.hpp"
#include "robotoc/dynamics/state_equation.hpp"
#include "robotoc/dynamics/contact_dynamics.hpp"
#include "robotoc/dynamics/switching_constraint.hpp"

#include <cassert>


namespace robotoc {

IntermediateStage::IntermediateStage(const std::shared_ptr<CostFunction>& cost,
                                     const std::shared_ptr<Constraints>& constraints,
                                     const std::shared_ptr<ContactSequence>& contact_sequence)
  : cost_(cost), 
    constraints_(constraints),
    contact_sequence_(contact_sequence) {
}


IntermediateStage::IntermediateStage()
  : cost_(), 
    constraints_(),
    contact_sequence_() {
}


OCPData IntermediateStage::createData(const Robot& robot) const {
  OCPData data;
  data.performance_index = PerformanceIndex();
  data.cost_data = cost_->createCostFunctionData(robot);
  data.constraints_data = constraints_->createConstraintsData(robot);
  data.state_equation_data = StateEquationData(robot);
  data.contact_dynamics_data = ContactDynamicsData(robot);
  data.switching_constraint_data = SwitchingConstraintData(robot);
  return data;
}


bool IntermediateStage::isFeasible(Robot& robot, const GridInfo& grid_info, 
                                   const SplitSolution& s, OCPData& data) const {
  const auto& contact_status = contact_sequence_->contactStatus(grid_info.contact_phase);
  return constraints_->isFeasible(robot, contact_status, data.constraints_data, s);
}


void IntermediateStage::initConstraints(Robot& robot, const GridInfo& grid_info, 
                                        const SplitSolution& s, OCPData& data) const {
  data.constraints_data.setTimeStage(grid_info.time_stage);
  const auto& contact_status = contact_sequence_->contactStatus(grid_info.contact_phase);
  constraints_->setSlackAndDual(robot, contact_status, data.constraints_data, s);
}


void IntermediateStage::evalOCP(Robot& robot, const GridInfo& grid_info, 
                                const SplitSolution& s, const SplitSolution& s_next, 
                                OCPData& data, SplitKKTResidual& kkt_residual) const {
  // setup computation
  const auto& contact_status = contact_sequence_->contactStatus(grid_info.contact_phase);
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_residual.setContactDimension(contact_status.dimf());
  kkt_residual.setZero();
  // eval cost and constraints
  data.performance_index.cost 
      = cost_->evalStageCost(robot, contact_status, data.cost_data, grid_info, s);
  constraints_->evalConstraint(robot, contact_status, data.constraints_data, s);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  // eval dynamics
  evalStateEquation(robot, grid_info.dt, s, s_next.q, s_next.v, kkt_residual);
  evalContactDynamics(robot, contact_status, s, data.contact_dynamics_data);
  if (grid_info.switching_constraint) {
    const auto& impulse_status = contact_sequence_->impulseStatus(grid_info.impulse_index+1);
    evalSwitchingConstraint(robot, impulse_status, data.switching_constraint_data,
                            grid_info.dt, grid_info.dt_next, s, kkt_residual);
  }
  // summarize evaluations
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
}


void IntermediateStage::evalKKT(Robot& robot, const GridInfo& grid_info, 
                                const Eigen::VectorXd& q_prev, 
                                const SplitSolution& s, const SplitSolution& s_next, 
                                OCPData& data, SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual) const {
  assert(q_prev.size() == robot.dimq());
  // setup computation
  const auto& contact_status = contact_sequence_->contactStatus(grid_info.contact_phase);
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactDimension(contact_status.dimf());
  kkt_residual.setContactDimension(contact_status.dimf());
  kkt_matrix.setZero();
  kkt_residual.setZero();
  // eval cost and constraints
  data.performance_index.cost 
      = cost_->quadratizeStageCost(robot, contact_status, data.cost_data,  
                                   grid_info, s, kkt_residual, kkt_matrix);
  kkt_residual.h  = (1.0/grid_info.dt) * data.performance_index.cost;
  kkt_matrix.hx   = (1.0/grid_info.dt) * kkt_residual.lx;
  kkt_matrix.hu   = (1.0/grid_info.dt) * kkt_residual.lu;
  kkt_matrix.ha   = (1.0/grid_info.dt) * kkt_residual.la;
  kkt_matrix.hf() = (1.0/grid_info.dt) * kkt_residual.lf();
  constraints_->linearizeConstraints(robot, contact_status, data.constraints_data, 
                                     s, kkt_residual);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  // eval dynamics
  linearizeStateEquation(robot, grid_info.dt, q_prev, s, s_next, 
                         data.state_equation_data, kkt_matrix, kkt_residual);
  linearizeContactDynamics(robot, contact_status, s, 
                           data.contact_dynamics_data, kkt_residual);
  if (grid_info.switching_constraint) {
    const auto& impulse_status = contact_sequence_->impulseStatus(grid_info.impulse_index+1);
    kkt_matrix.setSwitchingConstraintDimension(impulse_status.dimf());
    kkt_residual.setSwitchingConstraintDimension(impulse_status.dimf());
    linearizeSwitchingConstraint(robot, impulse_status, data.switching_constraint_data,
                                 grid_info.dt, grid_info.dt_next, s, 
                                 kkt_matrix, kkt_residual);
  }
  else {
    kkt_matrix.setSwitchingConstraintDimension(0);
    kkt_residual.setSwitchingConstraintDimension(0);
  }
  // summarize evaluations
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data.performance_index.dual_feasibility 
      = data.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data.performance_index.kkt_error = data.KKTError() + kkt_residual.KKTError();
  // Forms linear system
  constraints_->condenseSlackAndDual(contact_status, data.constraints_data, 
                                     kkt_matrix, kkt_residual);
  condenseContactDynamics(robot, contact_status, grid_info.dt, 
                          data.contact_dynamics_data, kkt_matrix, kkt_residual);
  correctLinearizeStateEquation(robot, grid_info.dt, s, s_next, 
                                data.state_equation_data, kkt_matrix, kkt_residual);
  kkt_residual.h        *= (1.0 / grid_info.N_phase);
  kkt_matrix.hx.array() *= (1.0 / grid_info.N_phase);
  kkt_matrix.hu.array() *= (1.0 / grid_info.N_phase);
  kkt_matrix.fx.array() *= (1.0 / grid_info.N_phase);
  kkt_matrix.Qtt        *= 1.0 / (grid_info.N_phase * grid_info.N_phase);
  kkt_matrix.Qtt_prev    = - kkt_matrix.Qtt;
}


void IntermediateStage::expandPrimal(const GridInfo& grid_info, OCPData& data, 
                                     SplitDirection& d) const {
  const auto& contact_status = contact_sequence_->contactStatus(grid_info.contact_phase);
  d.setContactDimension(contact_status.dimf());
  expandContactDynamicsPrimal(data.contact_dynamics_data, d);
  constraints_->expandSlackAndDual(contact_status, data.constraints_data, d);
}


void IntermediateStage::expandDual(const GridInfo& grid_info, OCPData& data,
                                   const SplitDirection& d_next, 
                                   SplitDirection& d, const double dts) const {
  assert(grid_info.dt > 0);
  expandContactDynamicsDual(grid_info.dt, dts, data.contact_dynamics_data, 
                            d_next, d);
  correctCostateDirection(data.state_equation_data, d);
}


double IntermediateStage::maxPrimalStepSize(const OCPData& data) const {
  return constraints_->maxSlackStepSize(data.constraints_data);
}


double IntermediateStage::maxDualStepSize(const OCPData& data) const {
  return constraints_->maxDualStepSize(data.constraints_data);
}


void IntermediateStage::updatePrimal(const Robot& robot, 
                                     const double primal_step_size, 
                                     const SplitDirection& d, 
                                     SplitSolution& s, OCPData& data) const {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(data.constraints_data, primal_step_size);
}


void IntermediateStage::updateDual(const double dual_step_size, 
                                   OCPData& data) const {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(data.constraints_data, dual_step_size);
}

} // namespace robotoc