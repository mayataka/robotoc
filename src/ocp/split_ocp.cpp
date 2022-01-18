#include "robotoc/ocp/split_ocp.hpp"

#include <cassert>


namespace robotoc {

SplitOCP::SplitOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                   const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    state_equation_(robot),
    contact_dynamics_(robot),
    switching_constraint_(robot),
    stage_cost_(0) {
}


SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    state_equation_(),
    contact_dynamics_(),
    switching_constraint_(),
    stage_cost_(0) {
}


SplitOCP::~SplitOCP() {
}


bool SplitOCP::isFeasible(Robot& robot, const ContactStatus& contact_status, 
                          const SplitSolution& s) {
  return constraints_->isFeasible(robot, contact_status, constraints_data_, s);
}


void SplitOCP::initConstraints(Robot& robot, const ContactStatus& contact_status, 
                               const int time_step, const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, contact_status, constraints_data_, s);
}


void SplitOCP::initConstraints(const SplitOCP& other) { 
  constraints_data_.copySlackAndDual(other.constraintsData());
}


const ConstraintsData& SplitOCP::constraintsData() const {
  return constraints_data_;
}


void SplitOCP::evalOCP(Robot& robot, const ContactStatus& contact_status,
                       const GridInfo& grid_info, const SplitSolution& s, 
                       const Eigen::VectorXd& q_next, 
                       const Eigen::VectorXd& v_next,
                       SplitKKTResidual& kkt_residual) {
  assert(q_next.size() == robot.dimq());
  assert(v_next.size() == robot.dimv());
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setZero();
  stage_cost_ = cost_->evalStageCost(robot, contact_status, cost_data_, 
                                     grid_info, s);
  constraints_->evalConstraint(robot, contact_status, constraints_data_, s);
  stage_cost_ += constraints_data_.logBarrier();
  state_equation_.evalStateEquation(robot, grid_info.dt, s, q_next, v_next, 
                                    kkt_residual);
  contact_dynamics_.evalContactDynamics(robot, contact_status, s);
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
  switching_constraint_.evalSwitchingConstraint(robot, impulse_status, 
                                                grid_info.dt, grid_info_next.dt, 
                                                s, sc_residual);
}


void SplitOCP::computeKKTResidual(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const GridInfo& grid_info, 
                                  const Eigen::VectorXd& q_prev, 
                                  const SplitSolution& s,
                                  const SplitSolution& s_next, 
                                  SplitKKTMatrix& kkt_matrix,
                                  SplitKKTResidual& kkt_residual) {
  computeKKTResidual_impl(robot, contact_status, grid_info, 
                          q_prev, s, s_next, kkt_matrix, kkt_residual);
}


void SplitOCP::computeKKTResidual(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const GridInfo& grid_info, 
                                  const Eigen::VectorXd& q_prev, 
                                  const SplitSolution& s,
                                  const ImpulseSplitSolution& s_next, 
                                  SplitKKTMatrix& kkt_matrix,
                                  SplitKKTResidual& kkt_residual) {
  computeKKTResidual_impl(robot, contact_status, grid_info, 
                          q_prev, s, s_next, kkt_matrix, kkt_residual);
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
  computeKKTResidual_impl(robot, contact_status, grid_info, 
                          q_prev, s, s_next, kkt_matrix, kkt_residual);
  switching_constraint_.linearizeSwitchingConstraint(robot, impulse_status, 
                                                     grid_info.dt, 
                                                     grid_info_next.dt, s, 
                                                     kkt_matrix, kkt_residual, 
                                                     sc_jacobian, sc_residual);
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
  computeKKTSystem_impl(robot, contact_status, grid_info, 
                        q_prev, s, s_next, kkt_matrix, kkt_residual);
}


void SplitOCP::computeKKTSystem(Robot& robot, 
                                const ContactStatus& contact_status,  
                                const GridInfo& grid_info, 
                                const Eigen::VectorXd& q_prev, 
                                const SplitSolution& s, 
                                const ImpulseSplitSolution& s_next,
                                SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual) {
  computeKKTSystem_impl(robot, contact_status, grid_info, 
                        q_prev, s, s_next, kkt_matrix, kkt_residual);
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
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeStageCost(robot, contact_status, cost_data_,  
                                           grid_info, s, kkt_residual, kkt_matrix);
  kkt_residual.h = (1.0/grid_info.dt) * stage_cost_;
  setHamiltonianDerivatives(grid_info.dt, kkt_matrix, kkt_residual);
  constraints_->linearizeConstraints(robot, contact_status, constraints_data_, 
                                     s, kkt_residual);
  stage_cost_ += constraints_data_.logBarrier();
  state_equation_.linearizeStateEquation(robot, grid_info.dt, q_prev, s, s_next, 
                                         kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, s,
                                             kkt_residual);
  switching_constraint_.linearizeSwitchingConstraint(robot, impulse_status, 
                                                     grid_info.dt, grid_info_next.dt, 
                                                     s, kkt_matrix, kkt_residual, 
                                                     sc_jacobian, sc_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual, sc_residual);
  constraints_->condenseSlackAndDual(contact_status, constraints_data_, 
                                     kkt_matrix, kkt_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, grid_info.dt, 
                                            kkt_matrix, kkt_residual);
  contact_dynamics_.condenseSwitchingConstraint(sc_jacobian, sc_residual, 
                                                kkt_matrix);
  state_equation_.correctLinearizedStateEquation(robot, grid_info.dt, s, s_next, 
                                                 kkt_matrix, kkt_residual);
}


void SplitOCP::computeInitialStateDirection(const Robot& robot, 
                                            const Eigen::VectorXd& q0, 
                                            const Eigen::VectorXd& v0, 
                                            const SplitSolution& s0, 
                                            SplitDirection& d0) const {
  state_equation_.computeInitialStateDirection(robot, q0, v0, s0, d0);
}


void SplitOCP::expandPrimal(const ContactStatus& contact_status, 
                            SplitDirection& d) {
  d.setContactStatus(contact_status);
  contact_dynamics_.expandPrimal(d);
  constraints_->expandSlackAndDual(contact_status, constraints_data_, d);
}


void SplitOCP::expandDual(const GridInfo& grid_info, 
                          const SplitDirection& d_next, 
                          SplitDirection& d, const double dts) {
  expandDual_impl(grid_info, d_next, d, dts);
}


void SplitOCP::expandDual(const GridInfo& grid_info, 
                          const ImpulseSplitDirection& d_next, 
                          SplitDirection& d, const double dts) {
  expandDual_impl(grid_info, d_next, d, dts);
}


void SplitOCP::expandDual(const GridInfo& grid_info, 
                          const SplitDirection& d_next, 
                          const SwitchingConstraintJacobian& sc_jacobian,
                          SplitDirection& d, const double dts) {
  assert(grid_info.dt > 0);
  contact_dynamics_.expandDual(grid_info.dt, dts, d_next, sc_jacobian, d);
  state_equation_.correctCostateDirection(d);
}


double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


void SplitOCP::updatePrimal(const Robot& robot, const double primal_step_size, 
                            const SplitDirection& d, SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


void SplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


double SplitOCP::KKTError(const SplitKKTResidual& kkt_residual) const {
  double err = 0;
  err += kkt_residual.KKTError();
  err += contact_dynamics_.KKTError();
  err += constraints_data_.KKTError();
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


double SplitOCP::stageCost() const {
  return stage_cost_;
} 


double SplitOCP::constraintViolation(
    const SplitKKTResidual& kkt_residual) const {
  double vio = 0;
  vio += kkt_residual.constraintViolation();
  vio += contact_dynamics_.constraintViolation();
  vio += constraints_data_.constraintViolation();
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