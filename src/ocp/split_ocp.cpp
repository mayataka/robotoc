#include "robotoc/ocp/split_ocp.hpp"
#include "robotoc/dynamics/state_equation.hpp"
#include "robotoc/dynamics/contact_dynamics.hpp"
#include "robotoc/dynamics/switching_constraint.hpp"

#include <cassert>


namespace robotoc {

SplitOCP::SplitOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints) {
  ocp_.cost = cost;
  ocp_.constraints = constraints;
  data_.cost_data = cost->createCostFunctionData(robot);
  data_.constraints_data = constraints->createConstraintsData(robot, 0);
  data_.state_equation_data = StateEquationData(robot);
  data_.contact_dynamics_data = ContactDynamicsData(robot);
  data_.switching_constraint_data = SwitchingConstraintData(robot);
}


SplitOCP::SplitOCP() 
  : ocp_(),
    data_() {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(Robot& robot, const ContactStatus& contact_status, 
                          const SplitSolution& s) {
  return ocp_.constraints->isFeasible(robot, contact_status, data_.constraints_data, s);
}


void SplitOCP::initConstraints(Robot& robot, const ContactStatus& contact_status, 
                               const int time_step, const SplitSolution& s) { 
  assert(time_step >= 0);
  data_.constraints_data.setTimeStage(time_step);
  ocp_.constraints->setSlackAndDual(robot, contact_status, data_.constraints_data, s);
}


void SplitOCP::initConstraints(const SplitOCP& other) { 
  data_.constraints_data.copySlackAndDual(other.constraintsData());
}


const ConstraintsData& SplitOCP::constraintsData() const {
  return data_.constraints_data;
}


void SplitOCP::evalOCP(Robot& robot, const ContactStatus& contact_status,
                       const GridInfo& grid_info, const SplitSolution& s, 
                       const Eigen::VectorXd& q_next, 
                       const Eigen::VectorXd& v_next,
                       SplitKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_residual.setContactDimension(contact_status.dimf());
  kkt_residual.setZero();
  data_.performance_index.cost 
      = ocp_.cost->evalStageCost(robot, contact_status, data_.cost_data, grid_info, s);
  ocp_.constraints->evalConstraint(robot, contact_status, data_.constraints_data, s);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  evalStateEquation(robot, grid_info.dt, s, q_next, v_next, kkt_residual);
  evalContactDynamics(robot, contact_status, s, data_.contact_dynamics_data);
}


void SplitOCP::evalOCP(Robot& robot, const ContactStatus& contact_status,
                       const GridInfo& grid_info, const SplitSolution& s, 
                       const Eigen::VectorXd& q_next, 
                       const Eigen::VectorXd& v_next,
                       SplitKKTResidual& kkt_residual,
                       const ImpulseStatus& impulse_status,
                       const GridInfo& grid_info_next, 
                       SwitchingConstraintResidual& sc_residual) {
  evalOCP(robot, contact_status, grid_info, s, q_next, v_next, kkt_residual);
  evalSwitchingConstraint(robot, impulse_status, data_.switching_constraint_data,
                          grid_info.dt, grid_info_next.dt, s, sc_residual);
}


void SplitOCP::computeKKTResidual(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const GridInfo& grid_info, 
                                  const Eigen::VectorXd& q_prev, 
                                  const SplitSolution& s,
                                  const SplitSolution& s_next, 
                                  SplitKKTMatrix& kkt_matrix,
                                  SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactDimension(contact_status.dimf());
  kkt_residual.setContactDimension(contact_status.dimf());
  kkt_residual.setZero();
  data_.performance_index.cost 
      = ocp_.cost->linearizeStageCost(robot, contact_status, data_.cost_data, 
                                  grid_info, s, kkt_residual);
  kkt_residual.h = (1.0/grid_info.dt) * data_.performance_index.cost;
  ocp_.constraints->linearizeConstraints(robot, contact_status, data_.constraints_data, 
                                     s, kkt_residual);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  linearizeStateEquation(robot, grid_info.dt, q_prev, s, s_next, 
                         data_.state_equation_data, kkt_matrix, kkt_residual);
  linearizeContactDynamics(robot, contact_status, s,
                           data_.contact_dynamics_data, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
}


void SplitOCP::computeKKTResidual(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const GridInfo& grid_info, 
                                  const Eigen::VectorXd& q_prev, 
                                  const SplitSolution& s,
                                  const SplitSolution& s_next, 
                                  SplitKKTMatrix& kkt_matrix,
                                  SplitKKTResidual& kkt_residual,
                                  const ImpulseStatus& impulse_status, 
                                  const GridInfo& grid_info_next, 
                                  SwitchingConstraintJacobian& sc_jacobian,
                                  SwitchingConstraintResidual& sc_residual) {
  computeKKTResidual(robot, contact_status, grid_info, 
                     q_prev, s, s_next, kkt_matrix, kkt_residual);
  linearizeSwitchingConstraint(robot, impulse_status, data_.switching_constraint_data,
                               grid_info.dt, grid_info_next.dt, s, 
                               kkt_matrix, kkt_residual, sc_jacobian, sc_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual, sc_residual);
}


void SplitOCP::computeKKTSystem(Robot& robot, 
                                const ContactStatus& contact_status,  
                                const GridInfo& grid_info, 
                                const Eigen::VectorXd& q_prev, 
                                const SplitSolution& s, 
                                const SplitSolution& s_next,
                                SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactDimension(contact_status.dimf());
  kkt_residual.setContactDimension(contact_status.dimf());
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data_.performance_index.cost = ocp_.cost->quadratizeStageCost(robot, contact_status, data_.cost_data,  
                                           grid_info, s, kkt_residual, kkt_matrix);
  kkt_residual.h = (1.0/grid_info.dt) * data_.performance_index.cost;
  setHamiltonianDerivatives(grid_info.dt, kkt_matrix, kkt_residual);
  ocp_.constraints->linearizeConstraints(robot, contact_status, data_.constraints_data, 
                                     s, kkt_residual);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  linearizeStateEquation(robot, grid_info.dt, q_prev, s, s_next, 
                         data_.state_equation_data, kkt_matrix, kkt_residual);
  linearizeContactDynamics(robot, contact_status, s, 
                           data_.contact_dynamics_data, kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
  ocp_.constraints->condenseSlackAndDual(contact_status, data_.constraints_data, 
                                     kkt_matrix, kkt_residual);
  condenseContactDynamics(robot, contact_status, grid_info.dt, 
                          data_.contact_dynamics_data, 
                          kkt_matrix, kkt_residual);
  correctLinearizeStateEquation(robot, grid_info.dt, s, s_next, 
                                data_.state_equation_data, kkt_matrix, kkt_residual);
}


void SplitOCP::computeKKTSystem(Robot& robot, 
                                const ContactStatus& contact_status, 
                                const GridInfo& grid_info, 
                                const Eigen::VectorXd& q_prev, 
                                const SplitSolution& s, 
                                const SplitSolution& s_next, 
                                SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual, 
                                const ImpulseStatus& impulse_status,
                                const GridInfo& grid_info_next, 
                                SwitchingConstraintJacobian& sc_jacobian,
                                SwitchingConstraintResidual& sc_residual) {
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactDimension(contact_status.dimf());
  kkt_residual.setContactDimension(contact_status.dimf());
  kkt_matrix.setZero();
  kkt_residual.setZero();
  data_.performance_index.cost = ocp_.cost->quadratizeStageCost(robot, contact_status, data_.cost_data,  
                                           grid_info, s, kkt_residual, kkt_matrix);
  kkt_residual.h = (1.0/grid_info.dt) * data_.performance_index.cost;
  setHamiltonianDerivatives(grid_info.dt, kkt_matrix, kkt_residual);
  ocp_.constraints->linearizeConstraints(robot, contact_status, data_.constraints_data, 
                                     s, kkt_residual);
  data_.performance_index.cost_barrier = data_.constraints_data.logBarrier();
  linearizeStateEquation(robot, grid_info.dt, q_prev, s, s_next, 
                         data_.state_equation_data, kkt_matrix, kkt_residual);
  linearizeContactDynamics(robot, contact_status, s, 
                           data_.contact_dynamics_data, kkt_residual);
  linearizeSwitchingConstraint(robot, impulse_status, data_.switching_constraint_data,
                               grid_info.dt, grid_info_next.dt, s, 
                               kkt_matrix, kkt_residual, sc_jacobian, sc_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual, sc_residual);
  ocp_.constraints->condenseSlackAndDual(contact_status, data_.constraints_data, 
                                     kkt_matrix, kkt_residual);
  kkt_matrix.setSwitchingConstraintDimension(0);
  kkt_residual.setSwitchingConstraintDimension(0);
  condenseContactDynamics(robot, contact_status, grid_info.dt, 
                          data_.contact_dynamics_data, kkt_matrix, kkt_residual);
  condenseContactDynamics(data_.contact_dynamics_data, sc_jacobian, sc_residual);
  correctLinearizeStateEquation(robot, grid_info.dt, s, s_next, 
                                data_.state_equation_data, kkt_matrix, kkt_residual);
}


void SplitOCP::computeInitialStateDirection(const Robot& robot, 
                                            const Eigen::VectorXd& q0, 
                                            const Eigen::VectorXd& v0, 
                                            const SplitSolution& s0, 
                                            SplitDirection& d0) const {
  ::robotoc::computeInitialStateDirection(robot, q0, v0, s0, 
                                          data_.state_equation_data, d0);
}


void SplitOCP::expandPrimal(const ContactStatus& contact_status, 
                            SplitDirection& d) {
  d.setContactDimension(contact_status.dimf());
  expandContactDynamicsPrimal(data_.contact_dynamics_data, d);
  ocp_.constraints->expandSlackAndDual(contact_status, data_.constraints_data, d);
}


void SplitOCP::expandDual(const GridInfo& grid_info, 
                          const SplitDirection& d_next, 
                          SplitDirection& d, const double dts) {
  assert(grid_info.dt > 0);
  expandContactDynamicsDual(grid_info.dt, dts, data_.contact_dynamics_data, 
                            d_next, d);
  correctCostateDirection(data_.state_equation_data, d);
}


void SplitOCP::expandDual(const GridInfo& grid_info, 
                          const SplitDirection& d_next, 
                          const SwitchingConstraintJacobian& sc_jacobian,
                          SplitDirection& d, const double dts) {
  assert(grid_info.dt > 0);
  expandContactDynamicsDual(grid_info.dt, dts, data_.contact_dynamics_data, 
                            sc_jacobian, d_next, d);
  correctCostateDirection(data_.state_equation_data, d);
}


double SplitOCP::maxPrimalStepSize() {
  return ocp_.constraints->maxSlackStepSize(data_.constraints_data);
}


double SplitOCP::maxDualStepSize() {
  return ocp_.constraints->maxDualStepSize(data_.constraints_data);
}


void SplitOCP::updatePrimal(const Robot& robot, const double primal_step_size, 
                            const SplitDirection& d, SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  ocp_.constraints->updateSlack(data_.constraints_data, primal_step_size);
}


void SplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  ocp_.constraints->updateDual(data_.constraints_data, dual_step_size);
}


double SplitOCP::KKTError(const SplitKKTResidual& kkt_residual) const {
  double err = 0;
  err += kkt_residual.KKTError();
  err += data_.contact_dynamics_data.KKTError();
  err += data_.constraints_data.KKTError();
  return err;
}


double SplitOCP::KKTError(
    const SplitKKTResidual& kkt_residual, 
    const SwitchingConstraintResidual& sc_residual) const {
  double err = 0;
  err += KKTError(kkt_residual);
  err += sc_residual.KKTError();
  return err;
}


double SplitOCP::stageCost(const bool include_cost_barrier) const {
  if (include_cost_barrier) {
    return data_.performance_index.cost + data_.performance_index.cost_barrier; 
  }
  else {
    return data_.performance_index.cost;
  }
} 


double SplitOCP::constraintViolation(
    const SplitKKTResidual& kkt_residual) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += data_.contact_dynamics_data.constraintViolation();
  vio += data_.constraints_data.constraintViolation();
  return vio;
}


double SplitOCP::constraintViolation(
    const SplitKKTResidual& kkt_residual, 
    const SwitchingConstraintResidual& sc_residual) const {
  double vio = 0;
  vio += constraintViolation(kkt_residual);
  vio += sc_residual.constraintViolation();
  return vio;
}


void SplitOCP::setHamiltonianDerivatives(const double dt, 
                                         SplitKKTMatrix& kkt_matrix, 
                                         const SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  const double coeff = 1.0 / dt;
  kkt_matrix.hx = coeff * kkt_residual.lx;
  kkt_matrix.hu = coeff * kkt_residual.lu;
  kkt_matrix.ha = coeff * kkt_residual.la;
  kkt_matrix.hf() = coeff * kkt_residual.lf();
}


void SplitOCP::correctSTOSensitivities(SplitKKTMatrix& kkt_matrix,
                                       SplitKKTResidual& kkt_residual,
                                       const int N_phase) {
  assert(N_phase > 0);
  const double coeff = 1.0 / N_phase;
  kkt_residual.h        *= coeff;
  kkt_matrix.hx.array() *= coeff;
  kkt_matrix.hu.array() *= coeff;
  kkt_matrix.fx.array() *= coeff;
  const double coeff2 = 1.0 / (N_phase*N_phase);
  kkt_matrix.Qtt        *= coeff2;
  kkt_matrix.Qtt_prev    = - kkt_matrix.Qtt;
}


void SplitOCP::correctSTOSensitivities(SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual, 
                                       SwitchingConstraintJacobian& sc_jacobian, 
                                       const int N_phase) {
  assert(N_phase > 0);
  correctSTOSensitivities(kkt_matrix, kkt_residual, N_phase);
  const double coeff = 1.0 / N_phase;
  sc_jacobian.Phit().array() *= coeff;
}

} // namespace robotoc