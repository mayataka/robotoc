#include "idocp/impulse/impulse_split_parnmpc.hpp"

#include <cassert>


namespace idocp {

inline ImpulseSplitParNMPC::ImpulseSplitParNMPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, -1)),
    impulse_dynamics_(robot) {
}


inline ImpulseSplitParNMPC::ImpulseSplitParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    impulse_dynamics_() {
}


inline ImpulseSplitParNMPC::~ImpulseSplitParNMPC() {
}


inline bool ImpulseSplitParNMPC::isFeasible(Robot& robot, 
                                            const ImpulseSplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void ImpulseSplitParNMPC::initConstraints(
    Robot& robot, const ImpulseSplitSolution& s) {
  constraints_data_ = constraints_->createConstraintsData(robot, -1);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void ImpulseSplitParNMPC::linearizeOCP(
    Robot& robot, const ImpulseStatus& impulse_status, const double t,  
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v);
  // condensing the impulse dynamics
  kkt_residual.setZero();
  kkt_matrix.setZero();
  cost_->computeImpulseCostDerivatives(robot, cost_data_, t, s, kkt_residual);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s,
                                             kkt_matrix, kkt_residual);
  cost_->computeImpulseCostHessian(robot, cost_data_, t, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, s, kkt_matrix, 
                                     kkt_residual);
  impulse_dynamics_.condenseImpulseDynamics(robot, impulse_status, 
                                            kkt_matrix, kkt_residual);
}


inline void ImpulseSplitParNMPC::computeCondensedPrimalDirection(
    Robot& robot, const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitSolution& s, ImpulseSplitDirection& d) {
  d.setImpulseStatusByDimension(s.dimf());
  impulse_dynamics_.computeCondensedPrimalDirection(robot, kkt_matrix, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline void ImpulseSplitParNMPC::computeCondensedDualDirection(
    const Robot& robot, ImpulseSplitDirection& d) {
  impulse_dynamics_.computeCondensedDualDirection(robot, d);
}

 
inline double ImpulseSplitParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double ImpulseSplitParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void ImpulseSplitParNMPC::updatePrimal(
    const Robot& robot, const double primal_step_size, 
    const ImpulseSplitDirection& d, ImpulseSplitSolution& s,
    const bool is_position_constraint_valid) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d, is_position_constraint_valid);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void ImpulseSplitParNMPC::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline void ImpulseSplitParNMPC::computeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double t,  
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next,
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v);
  cost_->computeImpulseCostDerivatives(robot, cost_data_, t, s, kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s,
                                             kkt_matrix, kkt_residual);
}


inline double ImpulseSplitParNMPC::squaredNormKKTResidual(
    const ImpulseSplitKKTResidual& kkt_residual) const {
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  error += kkt_residual.ldv.squaredNorm();
  error += kkt_residual.lf().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  error += impulse_dynamics_.squaredNormImpulseDynamicsResidual(kkt_residual);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double ImpulseSplitParNMPC::stageCost(Robot& robot, const double t, 
                                             const ImpulseSplitSolution& s, 
                                             const double primal_step_size) {
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  robot.updateKinematics(s.q, s.v);
  double cost = 0;
  cost += cost_->computeImpulseCost(robot, cost_data_, t, s);
  if (primal_step_size > 0) {
    cost += constraints_->costSlackBarrier(constraints_data_, primal_step_size);
  }
  else {
    cost += constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double ImpulseSplitParNMPC::constraintViolation(
    Robot& robot, const ImpulseStatus& impulse_status, const double t, 
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev,
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::ComputeImpulseBackwardEulerResidual(robot, q_prev, v_prev, s, 
                                                     kkt_residual);
  impulse_dynamics_.computeImpulseDynamicsResidual(robot, impulse_status, s, 
                                                   kkt_residual);
  double violation = 0;
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual);
  violation += impulse_dynamics_.l1NormImpulseDynamicsResidual(kkt_residual);
  return violation;
}

} // namespace idocp