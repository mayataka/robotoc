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
    constraints_data_(constraints->createConstraintsData(robot)),
    contact_dynamics_(robot),
    has_floating_base_(robot.has_floating_base()),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.max_point_contacts() > 0) {
    use_kinematics_ = true;
  }
}


inline SplitOCP::SplitOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    contact_dynamics_(),
    has_floating_base_(false),
    use_kinematics_(false) {
}


inline SplitOCP::~SplitOCP() {
}


inline bool SplitOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitOCP::initConstraints(Robot& robot, const int time_step, 
                                      const double dtau, 
                                      const SplitSolution& s) { 
  assert(time_step >= 0);
  assert(dtau > 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, dtau, s);
}


template <typename SplitSolutionType>
inline void SplitOCP::linearizeOCP(Robot& robot, 
                                   const ContactStatus& contact_status,  
                                   const double t, const double dtau, 
                                   const Eigen::VectorXd& q_prev, 
                                   const SplitSolution& s, 
                                   const SplitSolutionType& s_next,
                                   SplitKKTMatrix& kkt_matrix, 
                                   SplitKKTResidual& kkt_residual) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_matrix.setZero();
  kkt_residual.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dtau, s, 
                                             kkt_residual);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix, kkt_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dtau, 
                                            kkt_matrix, kkt_residual);
}


inline void SplitOCP::computeCondensedPrimalDirection(
    Robot& robot, const double dtau, const SplitSolution& s, 
    SplitDirection& d) {
  contact_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, 
                                             s, d);
}


template <typename SplitDirectionType>
inline void SplitOCP::computeCondensedDualDirection(
    const Robot& robot, const double dtau, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, const SplitDirectionType& d_next, 
    SplitDirection& d) {
  assert(dtau > 0);
  contact_dynamics_.computeCondensedDualDirection(robot, dtau, kkt_matrix,
                                                  kkt_residual, 
                                                  d_next.dgmm(), d);
}


inline double SplitOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitOCP::updatePrimal(const Robot& robot, 
                                   const double primal_step_size, 
                                   const double dtau, const SplitDirection& d, 
                                   SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  assert(dtau > 0);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


template <typename SplitSolutionType>
inline void SplitOCP::computeKKTResidual(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double t, const double dtau, 
                                         const Eigen::VectorXd& q_prev, 
                                         const SplitSolution& s,
                                         const SplitSolutionType& s_next, 
                                         SplitKKTMatrix& kkt_matrix,
                                         SplitKKTResidual& kkt_residual) {
  assert(dtau > 0);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dtau, s, 
                                             kkt_residual);
}


inline double SplitOCP::squaredNormKKTResidual(
    const SplitKKTResidual& kkt_residual, const double dtau) const {
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  error += kkt_residual.la.squaredNorm();
  error += kkt_residual.lf().squaredNorm();
  if (has_floating_base_) {
    error += kkt_residual.lu_passive.squaredNorm();
  }
  error += kkt_residual.lu().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  error += contact_dynamics_.squaredNormContactDynamicsResidual(dtau);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitOCP::stageCost(Robot& robot, const double t,  
                                  const double dtau, const SplitSolution& s, 
                                  const double primal_step_size) {
  assert(dtau > 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, dtau, s);
  if (primal_step_size > 0) {
    cost += constraints_->costSlackBarrier(constraints_data_, primal_step_size);
  }
  else {
    cost += constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitOCP::constraintViolation(Robot& robot, 
                                            const ContactStatus& contact_status, 
                                            const double t, const double dtau, 
                                            const SplitSolution& s, 
                                            const Eigen::VectorXd& q_next, 
                                            const Eigen::VectorXd& v_next,
                                            SplitKKTResidual& kkt_residual) {
  kkt_residual.setContactStatus(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  stateequation::ComputeForwardEulerResidual(robot, dtau, s, q_next, v_next, 
                                             kkt_residual);
  contact_dynamics_.computeContactDynamicsResidual(robot, contact_status, dtau, s);
  double violation = 0;
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual);
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  violation += contact_dynamics_.l1NormContactDynamicsResidual(dtau);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_OCP_HXX_