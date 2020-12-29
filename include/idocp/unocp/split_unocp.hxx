#ifndef IDOCP_SPLIT_UNOCP_HXX_
#define IDOCP_SPLIT_UNOCP_HXX_

#include "idocp/unocp/split_unocp.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline SplitUnOCP::SplitUnOCP(const Robot& robot, 
                              const std::shared_ptr<CostFunction>& cost,
                              const std::shared_ptr<Constraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot, 0)),
    unconstrained_dynamics_(robot),
    use_kinematics_(false) {
  if (cost_->useKinematics() || constraints_->useKinematics() 
                             || robot.maxPointContacts() > 0) {
    use_kinematics_ = true;
  }
  try {
    if (robot.hasFloatingBase()) {
      throw std::logic_error("robot has floating base: robot should have no constraints!");
    }
    if (robot.maxPointContacts() > 0) {
      throw std::logic_error("robot can have contacts: robot should have no constraints!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline SplitUnOCP::SplitUnOCP() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    unconstrained_dynamics_(),
    use_kinematics_(false) {
}


inline SplitUnOCP::~SplitUnOCP() {
}


inline bool SplitUnOCP::isFeasible(Robot& robot, const SplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitUnOCP::initConstraints(Robot& robot, const int time_step, 
                                      const SplitSolution& s) { 
  assert(time_step >= 0);
  constraints_data_ = constraints_->createConstraintsData(robot, time_step);
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


template <typename MatrixType1, typename MatrixType2>
inline void SplitUnOCP::linearizeOCP(
    Robot& robot, const double t, const double dtau, 
    const Eigen::VectorXd& q_prev, const SplitSolution& s, 
    const SplitSolution& s_next, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual, const Eigen::MatrixBase<MatrixType1>& Qaq, 
    const Eigen::MatrixBase<MatrixType2>& Qav) {
  assert(dtau >= 0);
  assert(q_prev.size() == robot.dimq());
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
  unconstrained_dynamics_.linearizeUnconstrainedDynamics(robot, dtau, s, 
                                                         kkt_residual);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix, kkt_residual);
  unconstrained_dynamics_.condenseUnconstrainedDynamics(
      robot, dtau, kkt_matrix, kkt_residual, 
      const_cast<Eigen::MatrixBase<MatrixType1>&>(Qaq),
      const_cast<Eigen::MatrixBase<MatrixType1>&>(Qav));
}


inline void SplitUnOCP::computeCondensedDirection(
    Robot& robot, const double dtau, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, const SplitSolution& s, 
    SplitDirection& d) {
  unconstrained_dynamics_.computeCondensedDirection(dtau, kkt_matrix, 
                                                    kkt_residual, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline double SplitUnOCP::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitUnOCP::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitUnOCP::updatePrimal(const Robot& robot, 
                                   const double primal_step_size, 
                                   const SplitDirection& d, 
                                   SplitSolution& s) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


inline void SplitUnOCP::updateDual(const double dual_step_size) {
  assert(dual_step_size > 0);
  assert(dual_step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline void SplitUnOCP::computeKKTResidual(Robot& robot, const double t, 
                                           const double dtau, 
                                           const Eigen::VectorXd& q_prev, 
                                           const SplitSolution& s,
                                           const SplitSolution& s_next, 
                                           SplitKKTMatrix& kkt_matrix,
                                           SplitKKTResidual& kkt_residual) {
  assert(dtau >= 0);
  assert(q_prev.size() == robot.dimq());
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_residual.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix, kkt_residual);
  unconstrained_dynamics_.linearizeUnconstrainedDynamics(robot, dtau, s, 
                                                         kkt_residual);
}


inline double SplitUnOCP::squaredNormKKTResidual(
    const SplitKKTResidual& kkt_residual, const double dtau) const {
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  error += kkt_residual.la.squaredNorm();
  error += kkt_residual.lu().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  error += unconstrained_dynamics_.squaredNormUnconstrainedDynamicsResidual(dtau);
  error += dtau * dtau * constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitUnOCP::stageCost(Robot& robot, const double t,  
                                    const double dtau, const SplitSolution& s, 
                                    const double primal_step_size) {
  assert(dtau >= 0);
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  double cost = 0;
  cost += cost_->computeStageCost(robot, cost_data_, t, dtau, s);
  if (primal_step_size > 0) {
    cost += dtau * constraints_->costSlackBarrier(constraints_data_, 
                                                  primal_step_size);
  }
  else {
    cost += dtau * constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitUnOCP::constraintViolation(Robot& robot, 
                                              const double t, const double dtau, 
                                              const SplitSolution& s, 
                                              const Eigen::VectorXd& q_next, 
                                              const Eigen::VectorXd& v_next,
                                              SplitKKTResidual& kkt_residual) {
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::ComputeForwardEulerResidual(robot, dtau, s, q_next, v_next, 
                                             kkt_residual);
  unconstrained_dynamics_.computeUnconstrainedDynamicsResidual(robot, s);
  double violation = 0;
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual);
  violation += unconstrained_dynamics_.l1NormUnconstrainedDynamicsResidual(dtau);
  violation += dtau * constraints_->l1NormPrimalResidual(constraints_data_);
  return violation;
}

} // namespace idocp

#endif // IDOCP_SPLIT_UNOCP_HXX_