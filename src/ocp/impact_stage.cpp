#include "robotoc/ocp/impact_stage.hpp"
#include "robotoc/dynamics/impact_state_equation.hpp"
#include "robotoc/dynamics/impact_dynamics.hpp"

#include <cassert>


namespace robotoc {

ImpactStage::ImpactStage(const std::shared_ptr<CostFunction>& cost,
                         const std::shared_ptr<Constraints>& constraints,
                         const std::shared_ptr<ContactSequence>& contact_sequence)
  : cost_(cost), 
    constraints_(constraints),
    contact_sequence_(contact_sequence) {
}


ImpactStage::ImpactStage()
  : cost_(), 
    constraints_(),
    contact_sequence_() {
}


OCPData ImpactStage::createData(const Robot& robot) const {
  OCPData data;
  data.performance_index = PerformanceIndex();
  data.cost_data = cost_->createCostFunctionData(robot);
  data.constraints_data = constraints_->createConstraintsData(robot);
  data.state_equation_data = StateEquationData(robot);
  data.contact_dynamics_data = ContactDynamicsData(robot);
  data.switching_constraint_data = SwitchingConstraintData(robot);
  return data;
}


bool ImpactStage::isFeasible(Robot& robot, const GridInfo& grid_info, 
                             const SplitSolution& s, OCPData& data) const {
  assert(grid_info.type == GridType::Impact);
  const auto& impact_status = contact_sequence_->impactStatus(grid_info.impact_index);
  return constraints_->isFeasible(robot, impact_status, data.constraints_data, s);
}


void ImpactStage::initConstraints(Robot& robot, const GridInfo& grid_info, 
                                  const SplitSolution& s, OCPData& data) const {
  assert(grid_info.type == GridType::Impact);
  data.constraints_data.setTimeStage(-1);
  const auto& impact_status = contact_sequence_->impactStatus(grid_info.impact_index);
  constraints_->setSlackAndDual(robot, impact_status, data.constraints_data, s);
}


void ImpactStage::evalOCP(Robot& robot, const GridInfo& grid_info, 
                          const SplitSolution& s, const SplitSolution& s_next, 
                          OCPData& data, SplitKKTResidual& kkt_residual) const {
  assert(grid_info.type == GridType::Impact);
  // setup computation
  const auto& impact_status = contact_sequence_->impactStatus(grid_info.impact_index);
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_residual.setContactDimension(impact_status.dimf());
  kkt_residual.setSwitchingConstraintDimension(0);
  kkt_residual.setZero();
  // eval cost and constraints
  data.performance_index.cost 
      = cost_->evalImpactCost(robot, impact_status, data.cost_data, grid_info, s);
  constraints_->evalConstraint(robot, impact_status, data.constraints_data, s);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  // eval dynamics
  evalImpactStateEquation(robot, s, s_next, kkt_residual);
  evalImpactDynamics(robot, impact_status, s, data.contact_dynamics_data);
  // summarize evaluations
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
}


void ImpactStage::evalKKT(Robot& robot, const GridInfo& grid_info, 
                          const Eigen::VectorXd& q_prev, const SplitSolution& s, 
                          const SplitSolution& s_next, OCPData& data, 
                          SplitKKTMatrix& kkt_matrix, 
                          SplitKKTResidual& kkt_residual) const {
  assert(grid_info.type == GridType::Impact);
  assert(q_prev.size() == robot.dimq());
  // setup computation
  const auto& impact_status = contact_sequence_->impactStatus(grid_info.impact_index);
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_matrix.setContactDimension(impact_status.dimf());
  kkt_matrix.setSwitchingConstraintDimension(0);
  kkt_residual.setContactDimension(impact_status.dimf());
  kkt_residual.setSwitchingConstraintDimension(0);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  // eval cost and constraints
  data.performance_index.cost 
      = cost_->quadratizeImpactCost(robot, impact_status, data.cost_data,  
                                     grid_info, s, kkt_residual, kkt_matrix);
  constraints_->linearizeConstraints(robot, impact_status, data.constraints_data, 
                                     s, kkt_residual);
  data.performance_index.cost_barrier = data.constraints_data.logBarrier();
  // eval dynamics
  linearizeImpactStateEquation(robot, q_prev, s, s_next, data.state_equation_data, 
                                kkt_matrix, kkt_residual);
  linearizeImpactDynamics(robot, impact_status, s, data.contact_dynamics_data, 
                           kkt_residual);
  // summarize evaluations
  data.performance_index.primal_feasibility 
      = data.primalFeasibility<1>() + kkt_residual.primalFeasibility<1>();
  data.performance_index.dual_feasibility 
      = data.dualFeasibility<1>() + kkt_residual.dualFeasibility<1>();
  data.performance_index.kkt_error = data.KKTError() + kkt_residual.KKTError();
  // Forms linear system
  constraints_->condenseSlackAndDual(impact_status, data.constraints_data, 
                                     kkt_matrix, kkt_residual);
  condenseImpactDynamics(robot, impact_status, data.contact_dynamics_data, 
                          kkt_matrix, kkt_residual);
  correctLinearizeImpactStateEquation(robot, s, s_next, data.state_equation_data, 
                                       kkt_matrix, kkt_residual);
}


void ImpactStage::expandPrimal(const GridInfo& grid_info, OCPData& data, 
                               SplitDirection& d) const {
  assert(grid_info.type == GridType::Impact);
  const auto& impact_status = contact_sequence_->impactStatus(grid_info.impact_index);
  d.setContactDimension(impact_status.dimf());
  d.setSwitchingConstraintDimension(0);
  expandImpactDynamicsPrimal(data.contact_dynamics_data, d);
  constraints_->expandSlackAndDual(impact_status, data.constraints_data, d);
}


void ImpactStage::expandDual(const GridInfo& grid_info, OCPData& data,
                             const SplitDirection& d_next, 
                             SplitDirection& d) const {
  assert(grid_info.type == GridType::Impact);
  expandImpactDynamicsDual(data.contact_dynamics_data, d_next, d);
  correctCostateDirection(data.state_equation_data, d);
}


double ImpactStage::maxPrimalStepSize(const OCPData& data) const {
  return constraints_->maxSlackStepSize(data.constraints_data);
}


double ImpactStage::maxDualStepSize(const OCPData& data) const {
  return constraints_->maxDualStepSize(data.constraints_data);
}


void ImpactStage::updatePrimal(const Robot& robot, const double primal_step_size, 
                               const SplitDirection& d, SplitSolution& s, 
                               OCPData& data) const {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d, true);
  constraints_->updateSlack(data.constraints_data, primal_step_size);
}


void ImpactStage::updateDual(const double dual_step_size, OCPData& data) const {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(data.constraints_data, dual_step_size);
}

} // namespace robotoc