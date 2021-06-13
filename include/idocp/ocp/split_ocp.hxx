#ifndef IDOCP_SPLIT_OCP_HXX_
#define IDOCP_SPLIT_OCP_HXX_

#include "idocp/ocp/split_ocp.hpp"

#include <cassert>

namespace idocp {

inline SplitOCP::SplitOCP(const Robot& robot, 
                          const std::shared_ptr<CostFunction>& cost,
                          const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    state_equation_(robot),
    contact_dynamics_(robot),
    stage_cost_(0) {
}


inline SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    state_equation_(),
    contact_dynamics_(),
    stage_cost_(0) {
}


inline SplitOCP::~SplitOCP() {
}


inline bool SplitOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitOCP::initConstraints(Robot& robot, const int time_step, 
                                      const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


template <typename SplitSolutionType>
inline void SplitOCP::computeKKTResidual(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double t, const double dt, 
                                         const Eigen::VectorXd& q_prev, 
                                         const SplitSolution& s,
                                         const SplitSolutionType& s_next, 
                                         SplitKKTMatrix& kkt_matrix,
                                         SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeStageCost(robot, cost_data_, t, dt, s, 
                                          kkt_residual);
  constraints_->linearizePrimalAndDualResidual(robot, constraints_data_, dt, s, 
                                               kkt_residual);
  state_equation_.linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
}


inline void SplitOCP::computeKKTResidual(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double t, const double dt, 
                                         const Eigen::VectorXd& q_prev, 
                                         const SplitSolution& s,
                                         const SplitSolution& s_next, 
                                         SplitKKTMatrix& kkt_matrix,
                                         SplitKKTResidual& kkt_residual,
                                         const ImpulseStatus& impulse_status, 
                                         const double dt_next,
                                         SplitSwitchingConstraintJacobian& switch_jacobian,
                                         SplitSwitchingConstraintResidual& switch_residual) {
  assert(dt > 0);
  assert(dt_next > 0);
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeStageCost(robot, cost_data_, t, dt, s, 
                                          kkt_residual);
  constraints_->linearizePrimalAndDualResidual(robot, constraints_data_, dt, s, 
                                               kkt_residual);
  state_equation_.linearizeForwardEuler(robot, dt, q_prev, s, s_next, 
                                        kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s, 
                                             kkt_residual);
  switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, dt, 
                                                    dt_next, s, kkt_residual, 
                                                    switch_jacobian, 
                                                    switch_residual);
}


template <typename SplitSolutionType>
inline void SplitOCP::computeKKTSystem(Robot& robot, 
                                       const ContactStatus& contact_status,  
                                       const double t, const double dt, 
                                       const Eigen::VectorXd& q_prev, 
                                       const SplitSolution& s, 
                                       const SplitSolutionType& s_next,
                                       SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeStageCost(robot, cost_data_, t, dt, s, 
                                           kkt_residual, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix, kkt_residual);
  state_equation_.linearizeForwardEulerLieDerivative(robot, dt, q_prev, s, s_next, 
                                                     kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s,
                                             kkt_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dt, 
                                            kkt_matrix, kkt_residual);
}


inline void SplitOCP::computeKKTSystem(Robot& robot, 
                                       const ContactStatus& contact_status, 
                                       const double t, const double dt, 
                                       const Eigen::VectorXd& q_prev, 
                                       const SplitSolution& s, 
                                       const SplitSolution& s_next, 
                                       SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual, 
                                       const ImpulseStatus& impulse_status,
                                       const double dt_next, 
                                       SplitSwitchingConstraintJacobian& switch_jacobian,
                                       SplitSwitchingConstraintResidual& switch_residual) {
  assert(dt > 0);
  assert(dt_next > 0);
  assert(q_prev.size() == robot.dimq());
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_matrix.setZero();
  kkt_residual.setZero();
  stage_cost_ = cost_->quadratizeStageCost(robot, cost_data_, t, dt, s, 
                                           kkt_residual, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dt, s, 
                                     kkt_matrix, kkt_residual);
  state_equation_.linearizeForwardEulerLieDerivative(robot, dt, q_prev, s, s_next, 
                                                     kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dt, s,
                                             kkt_residual);
  switchingconstraint::linearizeSwitchingConstraint(robot, impulse_status, dt, 
                                                    dt_next, s, kkt_residual, 
                                                    switch_jacobian,
                                                    switch_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dt, 
                                            kkt_matrix, kkt_residual);
  contact_dynamics_.condenseSwitchingConstraint(switch_jacobian, 
                                                switch_residual);
}


inline void SplitOCP::computeInitialStateDirection(const Robot& robot, 
                                                   const Eigen::VectorXd& q0, 
                                                   const Eigen::VectorXd& v0, 
                                                   const SplitSolution& s0, 
                                                   SplitDirection& d0) const {
  state_equation_.computeInitialStateDirection(robot, q0, v0, s0, d0);
}


inline void SplitOCP::expandPrimal(const SplitSolution& s, SplitDirection& d) {
  d.setContactStatusByDimension(s.dimf());
  contact_dynamics_.expandPrimal(d);
  constraints_->expandSlackAndDual(constraints_data_, s, d);
}


template <typename SplitDirectionType>
inline void SplitOCP::expandDual(const double dt, 
                                 const SplitDirectionType& d_next, 
                                 SplitDirection& d) {
  assert(dt > 0);
  contact_dynamics_.expandDual(dt, d_next, d);
  state_equation_.correctCostateDirection(d);
}


inline double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitOCP::updatePrimal(const Robot& robot, 
                                   const double primal_step_size, 
                                   const SplitDirection& d, 
                                   SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline double SplitOCP::squaredNormKKTResidual(
    const SplitKKTResidual& kkt_residual, const double dt) const {
  assert(dt > 0);
  double nrm = 0;
  nrm += kkt_residual.squaredNormKKTResidual();
  nrm += contact_dynamics_.squaredNormKKTResidual(dt);
  nrm += (dt*dt) * constraints_data_.squaredNormKKTResidual();
  return nrm;
}


inline double SplitOCP::stageCost(Robot& robot, const double t,  
                                  const double dt, const SplitSolution& s, 
                                  const double primal_step_size) {
  assert(dt > 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  robot.updateKinematics(s.q, s.v, s.a);
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


inline double SplitOCP::constraintViolation(Robot& robot, 
                                            const ContactStatus& contact_status, 
                                            const double t, const double dt, 
                                            const SplitSolution& s, 
                                            const Eigen::VectorXd& q_next, 
                                            const Eigen::VectorXd& v_next,
                                            SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  kkt_residual.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  state_equation_.computeForwardEulerResidual(robot, dt, s, q_next, v_next, 
                                              kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, s);
  double violation = 0;
  violation += kkt_residual.l1NormConstraintViolation();
  violation += dt * contact_dynamics_.l1NormConstraintViolation();
  violation += dt * constraints_data_.l1NormConstraintViolation();
  return violation;
}


inline double SplitOCP::constraintViolation(
    Robot& robot, const ContactStatus& contact_status, 
    const double t, const double dt, const SplitSolution& s, 
    const Eigen::VectorXd& q_next, const Eigen::VectorXd& v_next, 
    SplitKKTResidual& kkt_residual, const ImpulseStatus& impulse_status, 
    const double dt_next, SplitSwitchingConstraintResidual& switch_residual) {
  assert(dt > 0);
  assert(dt_next > 0);
  kkt_residual.setContactStatus(contact_status);
  robot.updateKinematics(s.q, s.v, s.a);
  state_equation_.computeForwardEulerResidual(robot, dt, s, q_next, v_next, 
                                              kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, s);
  switchingconstraint::computeSwitchingConstraintResidual(robot, impulse_status,  
                                                          dt, dt_next, s, 
                                                          switch_residual);
  double violation = 0;
  violation += kkt_residual.l1NormConstraintViolation();
  violation += dt * contact_dynamics_.l1NormConstraintViolation();
  violation += dt * constraints_data_.l1NormConstraintViolation();
  violation += switch_residual.l1NormConstraintViolation();
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_OCP_HXX_