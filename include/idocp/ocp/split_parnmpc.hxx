#ifndef IDOCP_SPLIT_PARNMPC_HXX_ 
#define IDOCP_SPLIT_PARNMPC_HXX_

#include "idocp/ocp/split_parnmpc.hpp"
#include "idocp/ocp/switching_constraint.hpp"

#include <cassert>


namespace idocp {

inline SplitParNMPC::SplitParNMPC(
    const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
    const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(robot),
    constraints_(constraints),
    constraints_data_(),
    contact_dynamics_(robot),
    has_floating_base_(robot.hasFloatingBase()),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.maxPointContacts() > 0) {
    use_kinematics_ = true;
  }
}


inline SplitParNMPC::SplitParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    contact_dynamics_(),
    has_floating_base_(false),
    use_kinematics_(false) {
}


inline SplitParNMPC::~SplitParNMPC() {
}


inline bool SplitParNMPC::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitParNMPC::initConstraints(Robot& robot, const int time_step, 
                                          const SplitSolution& s) {
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


template <typename SplitSolutionType>
inline void SplitParNMPC::linearizeOCP(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       const double t, const double dt, 
                                       const Eigen::VectorXd& q_prev, 
                                       const Eigen::VectorXd& v_prev, 
                                       const SplitSolution& s, 
                                       const SplitSolutionType& s_next, 
                                       SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setImpulseStatus();
  kkt_residual.setImpulseStatus();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, 
                                     kkt_residual);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s,
                                    kkt_residual);
  stateequation::linearizeBackwardEuler(robot, dt, q_prev, v_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  stateequation::condenseBackwardEuler(robot, dt, q_prev, s, 
                                       kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
  cost_->computeStageCostHessian(robot, cost_data_, t, dt, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix, kkt_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dt, 
                                            kkt_matrix, kkt_residual, false);
}


template <typename SplitSolutionType>
inline void SplitParNMPC::linearizeOCP(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       const double t, const double dt, 
                                       const Eigen::VectorXd& q_prev, 
                                       const Eigen::VectorXd& v_prev, 
                                       const SplitSolution& s, 
                                       const SplitSolutionType& s_next, 
                                       SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual,
                                       const ImpulseStatus& impulse_status) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, 
                                     kkt_residual);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s,
                                    kkt_residual);
  stateequation::linearizeBackwardEuler(robot, dt, q_prev, v_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  stateequation::condenseBackwardEuler(robot, dt, q_prev, s, 
                                       kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
  switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, s, 
                                                    kkt_matrix, kkt_residual);
  cost_->computeStageCostHessian(robot, cost_data_, t, dt, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix, kkt_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dt, 
                                            kkt_matrix, kkt_residual, false);
}


inline void SplitParNMPC::computeCondensedPrimalDirection(
    Robot& robot, const double dt, const SplitSolution& s, SplitDirection& d) {
  assert(dt > 0);
  d.setContactStatusByDimension(s.dimf());
  contact_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline void SplitParNMPC::computeCondensedDualDirection(
    const Robot& robot, const double dt, const SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, SplitDirection& d) {
  assert(dt > 0);
  contact_dynamics_.computeCondensedDualDirection(robot, dt, kkt_matrix,
                                                  kkt_residual, d.dgmm(), d);
  stateequation::correctCostateDirectionBackwardEuler(robot, kkt_matrix, 
                                                      kkt_residual, d.dlmd());
}


inline double SplitParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitParNMPC::updatePrimal(const Robot& robot, 
                                       const double primal_step_size, 
                                       const SplitDirection& d, 
                                       SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitParNMPC::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


template <typename SplitSolutionType>
inline void SplitParNMPC::computeKKTResidual(
    Robot& robot, const ContactStatus& contact_status, const double t, 
    const double dt, const Eigen::VectorXd& q_prev, 
    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
    const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setImpulseStatus();
  kkt_residual.setImpulseStatus();
  kkt_residual.setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, 
                                     kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s,
                                    kkt_residual);
  stateequation::linearizeBackwardEuler(robot, dt, q_prev, v_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
}


template <typename SplitSolutionType>
inline void SplitParNMPC::computeKKTResidual(
    Robot& robot, const ContactStatus& contact_status, const double t, 
    const double dt, const Eigen::VectorXd& q_prev, 
    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
    const SplitSolutionType& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, const ImpulseStatus& impulse_status) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dt, s, 
                                     kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dt, s,
                                    kkt_residual);
  stateequation::linearizeBackwardEuler(robot, dt, q_prev, v_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
  switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, s, 
                                                    kkt_matrix, kkt_residual);
}


inline double SplitParNMPC::squaredNormKKTResidual(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  assert(dt > 0);
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  error += kkt_residual.la.squaredNorm();
  error += kkt_residual.lf().squaredNorm();
  if (has_floating_base_) {
    error += kkt_residual.lu_passive.squaredNorm();
  }
  error += kkt_residual.lu().squaredNorm();
  error += stateequation::squaredNormStateEuqationResidual(kkt_residual);
  error += contact_dynamics_.squaredNormContactDynamicsResidual(dt);
  error += dt * dt * constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  error += switchingconstraint::squaredNormSwitchingConstraintResidual(kkt_residual);
  return error;
}


inline double SplitParNMPC::stageCost(Robot& robot, const double t, 
                                      const double dt, const SplitSolution& s,
                                      const double primal_step_size) {
  assert(dt > 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->computeStageCost(robot, cost_data_, t, dt, s);
  if (primal_step_size > 0) {
    cost += dt * constraints_->costSlackBarrier(constraints_data_, 
                                                  primal_step_size);
  }
  else {
    cost += dt * constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitParNMPC::constraintViolation(
    Robot& robot, const ContactStatus& contact_status, const double t, 
    const double dt, const Eigen::VectorXd& q_prev, 
    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setImpulseStatus();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::computeBackwardEulerResidual(robot, dt, q_prev, v_prev, s, 
                                              kkt_residual);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, s);
  double violation = 0;
  violation += stateequation::l1NormStateEuqationResidual(kkt_residual);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dt);
  violation += dt * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}


inline double SplitParNMPC::constraintViolation(
    Robot& robot, const ContactStatus& contact_status, 
    const double t, const double dt, const Eigen::VectorXd& q_prev, 
    const Eigen::VectorXd& v_prev, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual, const ImpulseStatus& impulse_status) {
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setImpulseStatus(impulse_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::computeBackwardEulerResidual(robot, dt, q_prev, v_prev, s, 
                                              kkt_residual);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, s);
  switchingconstraint::computeSwitchingConstraintResidual(robot, impulse_status, 
                                                          kkt_residual);
  double violation = 0;
  violation += stateequation::l1NormStateEuqationResidual(kkt_residual);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dt);
  violation += dt * constraints_->l1NormPrimalResidual(constraints_data_);
  violation += switchingconstraint::l1NormSwitchingConstraintResidual(kkt_residual);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_PARNMPC_HXX_