#ifndef IDOCP_SPLIT_OCP_HXX_
#define IDOCP_SPLIT_OCP_HXX_

#include "idocp/ocp/split_ocp.hpp"

#include <cassert>

namespace idocp {

template <typename SplitSolutionType>
inline void SplitOCP::linearizeOCP(Robot& robot, 
                                   const ContactStatus& contact_status,  
                                   const double t, const double dtau, 
                                   const Eigen::VectorXd& q_prev, 
                                   const SplitSolution& s, 
                                   const SplitSolutionType& s_next) {
  assert(dtau > 0);
  assert(q_prev.size() == robot.dimq());
  setContactStatusForKKT(contact_status);
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  kkt_residual_.setZero();
  kkt_matrix_.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dtau, s, 
                                             kkt_matrix_, kkt_residual_);
  cost_->computeStageCostHessian(robot, cost_data_, t, dtau, s, kkt_matrix_);
  constraints_->condenseSlackAndDual(robot, constraints_data_, dtau, s, 
                                     kkt_matrix_, kkt_residual_);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, dtau, 
                                            kkt_matrix_, kkt_residual_);
}


template <typename MatrixType1, typename MatrixType2>
inline void SplitOCP::backwardStateConstraintFactorization(
    const Eigen::MatrixBase<MatrixType1>& T_next, const double dtau,
    const Eigen::MatrixBase<MatrixType2>& T) const {
  riccati_factorizer_.backwardStateConstraintFactorization(
      kkt_matrix_, T_next, dtau, 
      const_cast<Eigen::MatrixBase<MatrixType2>&> (T));
}


template <typename VectorType>
inline void SplitOCP::computePrimalDirection(
    Robot& robot, const double dtau, const RiccatiFactorization& riccati, 
    const SplitSolution& s, const Eigen::MatrixBase<VectorType>& dx0, 
    SplitDirection& d, const bool exist_state_constraint) {
  riccati_factorizer_.computeStateDirection(riccati, dx0, d, 
                                            exist_state_constraint);
  riccati_factorizer_.computeCostateDirection(riccati, d, 
                                              exist_state_constraint);
  riccati_factorizer_.computeControlInputDirection(riccati, d, 
                                                   exist_state_constraint);
  contact_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, dtau, 
                                             s, d);
}


template <typename SplitDirectionType>
inline void SplitOCP::computeDualDirection(
    Robot& robot, const double dtau, const SplitDirectionType& d_next, 
    SplitDirection& d) {
  assert(dtau > 0);
  contact_dynamics_.computeCondensedDualDirection(robot, dtau, kkt_matrix_,
                                                  kkt_residual_, 
                                                  d_next.dgmm(), d);
}


template <typename SplitSolutionType>
inline void SplitOCP::computeKKTResidual(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double t, const double dtau, 
                                         const Eigen::VectorXd& q_prev, 
                                         const SplitSolution& s,
                                         const SplitSolutionType& s_next) {
  assert(dtau > 0);
  setContactStatusForKKT(contact_status);
  kkt_residual_.setZero();
  if (use_kinematics_) {
    robot.updateKinematics(s.q, s.v, s.a);
  }
  cost_->computeStageCostDerivatives(robot, cost_data_, t, dtau, s, 
                                     kkt_residual_);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  constraints_->augmentDualResidual(robot, constraints_data_, dtau, s,
                                    kkt_residual_);
  stateequation::LinearizeForwardEuler(robot, dtau, q_prev, s, s_next, 
                                       kkt_matrix_, kkt_residual_);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, dtau, s, 
                                             kkt_matrix_, kkt_residual_);
}

} // namespace idocp

#endif // IDOCP_SPLIT_OCP_HXX_