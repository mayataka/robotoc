#include "robotoc/ocp/impulse_split_ocp.hpp"
#include "robotoc/dynamics/impulse_state_equation.hpp"
#include "robotoc/dynamics/impulse_dynamics.hpp"

#include <cassert>

namespace robotoc {

ImpulseSplitOCP::ImpulseSplitOCP(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) {
  ocp_.cost = cost;
  ocp_.constraints = constraints;
  data_.cost_data = cost->createCostFunctionData(robot);
  data_.constraints_data = constraints->createConstraintsData(robot, 0);
  data_.state_equation_data = StateEquationData(robot);
  data_.contact_dynamics_data = ContactDynamicsData(robot);
  data_.switching_constraint_data = SwitchingConstraintData(robot);
}


ImpulseSplitOCP::ImpulseSplitOCP() 
  : ocp_(),
    data_() {
}


ImpulseSplitOCP::~ImpulseSplitOCP() {
}


bool ImpulseSplitOCP::isFeasible(Robot& robot, 
                                 const ImpulseStatus& impulse_status,
                                 const SplitSolution& s) {
  return ocp_.constraints->isFeasible(robot, impulse_status, data_.constraints_data, s);
}


void ImpulseSplitOCP::initConstraints(Robot& robot,
                                      const ImpulseStatus& impulse_status,
                                      const SplitSolution& s) { 
  data_.constraints_data.setTimeStage(-1);
  ocp_.constraints->setSlackAndDual(robot, impulse_status, data_.constraints_data, s);
}


void ImpulseSplitOCP::initConstraints(const ImpulseSplitOCP& other) { 
  data_.constraints_data.copySlackAndDual(other.constraintsData());
}


const ConstraintsData& ImpulseSplitOCP::constraintsData() const {
  return data_.constraints_data;
}


void ImpulseSplitOCP::evalOCP(Robot& robot, const ImpulseStatus& impulse_status, 
                              const GridInfo& grid_info, 
                              const SplitSolution& s, 
                              const Eigen::VectorXd& q_next, 
                              const Eigen::VectorXd& v_next, 
                              SplitKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  kkt_residual.setContactDimension(impulse_status.dimf());
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v+s.dv);
  data_.performance_index.cost = ocp_.cost->evalImpulseCost(robot, impulse_status, data_.cost_data, 
                                       grid_info, s);
  ocp_.constraints->evalConstraint(robot, impulse_status, data_.constraints_data, s);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  evalImpulseStateEquation(robot, s, q_next, v_next, kkt_residual);
  evalImpulseDynamics(robot, impulse_status, s, data_.contact_dynamics_data);
}


void ImpulseSplitOCP::computeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const GridInfo& grid_info, 
    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_matrix.setContactDimension(impulse_status.dimf());
  kkt_residual.setContactDimension(impulse_status.dimf());
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data_.performance_index.cost = ocp_.cost->linearizeImpulseCost(robot, impulse_status, data_.cost_data, 
                                            grid_info, s, kkt_residual);
  ocp_.constraints->linearizeConstraints(robot, impulse_status, data_.constraints_data, 
                                     s, kkt_residual);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  linearizeImpulseStateEquation(robot, q_prev, s, s_next, 
                                data_.state_equation_data, kkt_matrix, kkt_residual);
  linearizeImpulseDynamics(robot, impulse_status, s, data_.contact_dynamics_data, kkt_residual);
}


void ImpulseSplitOCP::computeKKTSystem(
    Robot& robot, const ImpulseStatus& impulse_status, const GridInfo& grid_info,
    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v+s.dv);
  kkt_matrix.setContactDimension(impulse_status.dimf());
  kkt_residual.setContactDimension(impulse_status.dimf());
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data_.performance_index.cost = ocp_.cost->quadratizeImpulseCost(robot, impulse_status, data_.cost_data,  
                                             grid_info, s, kkt_residual, kkt_matrix);
  ocp_.constraints->linearizeConstraints(robot, impulse_status, data_.constraints_data, 
                                     s, kkt_residual);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  linearizeImpulseStateEquation(robot, q_prev, s, s_next, 
                                data_.state_equation_data, kkt_matrix, kkt_residual);
  linearizeImpulseDynamics(robot, impulse_status, s, data_.contact_dynamics_data, kkt_residual);
  ocp_.constraints->condenseSlackAndDual(impulse_status, data_.constraints_data, 
                                     kkt_matrix, kkt_residual);
  condenseImpulseDynamics(robot, impulse_status, data_.contact_dynamics_data, kkt_matrix, kkt_residual);
  correctLinearizeImpulseStateEquation(robot, s, s_next, 
                                       data_.state_equation_data, kkt_matrix, kkt_residual);
}


void ImpulseSplitOCP::expandPrimal(const ImpulseStatus& impulse_status, 
                                   SplitDirection& d) {
  d.setContactDimension(impulse_status.dimf());
  expandImpulseDynamicsPrimal(data_.contact_dynamics_data, d);
  ocp_.constraints->expandSlackAndDual(impulse_status, data_.constraints_data, d);
}


void ImpulseSplitOCP::expandDual(const SplitDirection& d_next, 
                                 SplitDirection& d) {
  expandImpulseDynamicsDual(data_.contact_dynamics_data, d_next, d);
  correctCostateDirection(data_.state_equation_data, d);
}


double ImpulseSplitOCP::maxPrimalStepSize() {
  return ocp_.constraints->maxSlackStepSize(data_.constraints_data);
}


double ImpulseSplitOCP::maxDualStepSize() {
  return ocp_.constraints->maxDualStepSize(data_.constraints_data);
}


void ImpulseSplitOCP::updatePrimal(const Robot& robot, 
                                   const double primal_step_size, 
                                   const SplitDirection& d, 
                                   SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d, true);
  ocp_.constraints->updateSlack(data_.constraints_data, primal_step_size);
}


void ImpulseSplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  ocp_.constraints->updateDual(data_.constraints_data, dual_step_size);
}


double ImpulseSplitOCP::KKTError(
    const SplitKKTResidual& kkt_residual) const {
  double err = 0;
  err += kkt_residual.KKTError();
  err += data_.contact_dynamics_data.KKTError();
  err += data_.constraints_data.KKTError();
  return err;
}


double ImpulseSplitOCP::stageCost(const bool include_cost_barrier) const {
  if (include_cost_barrier) {
    return data_.performance_index.cost + data_.performance_index.cost_barrier; 
  }
  else {
    return data_.performance_index.cost;
  }
}


double ImpulseSplitOCP::constraintViolation(
    const SplitKKTResidual& kkt_residual) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += data_.constraints_data.constraintViolation();
  vio += data_.contact_dynamics_data.constraintViolation();
  return vio;
}

} // namespace robotoc