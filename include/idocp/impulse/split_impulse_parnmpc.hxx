#include "idocp/impulse/split_impulse_parnmpc.hpp"

#include <cassert>


namespace idocp {

inline SplitImpulseParNMPC::SplitImpulseParNMPC(
    const Robot& robot, const std::shared_ptr<ImpulseCostFunction>& cost, 
    const std::shared_ptr<ImpulseConstraints>& constraints) 
  : cost_(cost),
    cost_data_(cost->createCostFunctionData(robot)),
    constraints_(constraints),
    constraints_data_(constraints->createConstraintsData(robot)),
    impulse_dynamics_(robot) {
}


inline SplitImpulseParNMPC::SplitImpulseParNMPC() 
  : cost_(),
    cost_data_(),
    constraints_(),
    constraints_data_(),
    impulse_dynamics_() {
}


inline SplitImpulseParNMPC::~SplitImpulseParNMPC() {
}


inline bool SplitImpulseParNMPC::isFeasible(Robot& robot, 
                                            const ImpulseSplitSolution& s) {
  return constraints_->isFeasible(robot, constraints_data_, s);
}


inline void SplitImpulseParNMPC::initConstraints(
    Robot& robot, const ImpulseSplitSolution& s) {
  constraints_->setSlackAndDual(robot, constraints_data_, s);
}


inline void SplitImpulseParNMPC::linearizeOCP(
    Robot& robot, const ImpulseStatus& impulse_status, const double t,  
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual,
    const bool is_state_constraint_valid) {
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v);
  // condensing the impulse dynamics
  kkt_residual.setZero();
  kkt_matrix.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s,
                                             kkt_matrix, kkt_residual,
                                             is_state_constraint_valid);
  cost_->computeStageCostHessian(robot, cost_data_, t, s, kkt_matrix);
  constraints_->condenseSlackAndDual(robot, constraints_data_, s, kkt_matrix, 
                                     kkt_residual);
  impulse_dynamics_.condenseImpulseDynamics(robot, impulse_status, 
                                            kkt_matrix, kkt_residual);
}



// void SplitImpulseParNMPC::coarseUpdate(Robot& robot, 
//                                        const ContactStatus& contact_status,
//                                        const double t, 
//                                        const Eigen::VectorXd& q_prev, 
//                                        const Eigen::VectorXd& v_prev,
//                                        const ImpulseSplitSolution& s,
//                                        const SplitSolution& s_next,
//                                        const Eigen::MatrixXd& aux_mat_next_old,
//                                        ImpulseSplitDirection& d, 
//                                        ImpulseSplitSolution& s_new_coarse) {
//   assert(q_prev.size() == robot.dimq());
//   assert(v_prev.size() == robot.dimv());
//   assert(aux_mat_next_old.rows() == 2*robot.dimv());
//   assert(aux_mat_next_old.cols() == 2*robot.dimv());
//   setContactStatusForKKT(contact_status);
//   robot.updateKinematics(s.q, s.v);
//   // condensing the impulse dynamics
//   kkt_matrix_.setZero();
//   kkt_residual_.setZero();
//   cost_->computeStageCostHessian(robot, cost_data_, t, s, kkt_matrix_);
//   cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual_);
//   stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
//                                                kkt_matrix_, kkt_residual_);
//   impulse_dynamics_.condenseImpulseDynamics(robot, contact_status, s, 
//                                             kkt_matrix_, kkt_residual_);
//   // coarse update of the solution
//   kkt_matrix_.Qxx().noalias() += aux_mat_next_old;
//   kkt_matrix_.symmetrize(); 
//   dimKKT_ = kkt_residual_.dimKKT();
//   kkt_matrix_.invert(kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_));
//   d.split_direction() = kkt_matrix_inverse_.topLeftCorner(dimKKT_, dimKKT_)
//                           * kkt_residual_.KKT_residual();
//   s_new_coarse.lmd = s.lmd - d.dlmd();
//   s_new_coarse.gmm = s.gmm - d.dgmm();
//   s_new_coarse.mu_stack() = s.mu_stack() - d.dmu();
//   s_new_coarse.f_stack() = s.f_stack() - d.df();
//   robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse.q);
//   s_new_coarse.v = s.v - d.dv();
// }


// void SplitImpulseParNMPC::getAuxiliaryMatrix(Eigen::MatrixXd& aux_mat) const {
//   assert(aux_mat.rows() == dimx_);
//   assert(aux_mat.cols() == dimx_);
//   aux_mat = - kkt_matrix_inverse_.topLeftCorner(dimx_, dimx_);
// }


// void SplitImpulseParNMPC::backwardCorrectionSerial(
//     const Robot& robot, const SplitSolution& s_next, 
//     const SplitSolution& s_next_new, ImpulseSplitSolution& s_new) {
//   x_res_.head(robot.dimv()) = s_next_new.lmd - s_next.lmd;
//   x_res_.tail(robot.dimv()) = s_next_new.gmm - s_next.gmm;
//   dx_.noalias() 
//       = kkt_matrix_inverse_.block(0, dimKKT_-dimx_, dimx_, dimx_) * x_res_;
//   s_new.lmd.noalias() -= dx_.head(robot.dimv());
//   s_new.gmm.noalias() -= dx_.tail(robot.dimv());
// }


// void SplitImpulseParNMPC::backwardCorrectionParallel(
//     const Robot& robot, ImpulseSplitDirection& d, 
//     ImpulseSplitSolution& s_new) const {
//   d.split_direction().tail(dimKKT_-dimx_).noalias()
//       = kkt_matrix_inverse_.block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) 
//           * x_res_;
//   s_new.mu_stack().noalias() -= d.dmu();
//   s_new.f_stack().noalias() -= d.df();
//   robot.integrateConfiguration(d.dq(), -1, s_new.q);
//   s_new.v.noalias() -= d.dv();
// }


// void SplitImpulseParNMPC::forwardCorrectionSerial(
//     const Robot& robot, const SplitSolution& s_prev, 
//     const SplitSolution& s_prev_prev, ImpulseSplitSolution& s_new) {
//   robot.subtractConfiguration(s_prev_prev.q, s_prev.q, 
//                               x_res_.head(robot.dimv()));
//   x_res_.tail(robot.dimv()) = s_prev_prev.v - s_prev.v;
//   dx_.noalias() = kkt_matrix_inverse_.block(dimKKT_-dimx_, 0, dimx_, dimx_) 
//                     * x_res_;
//   robot.integrateConfiguration(dx_.head(robot.dimv()), -1, s_new.q);
//   s_new.v.noalias() -= dx_.tail(robot.dimv());
// }


// void SplitImpulseParNMPC::forwardCorrectionParallel(
//     const Robot& robot, ImpulseSplitDirection& d, 
//     ImpulseSplitSolution& s_new) const {
//   d.split_direction().head(dimKKT_-dimx_).noalias()
//       = kkt_matrix_inverse_.topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
//   s_new.lmd.noalias() -= d.dlmd();
//   s_new.gmm.noalias() -= d.dgmm();
//   s_new.mu_stack().noalias() -= d.dmu();
//   s_new.f_stack().noalias() -= d.df();
// }


// void SplitImpulseParNMPC::computePrimalAndDualDirection(
//     Robot& robot, const ImpulseSplitSolution& s, 
//     const ImpulseSplitSolution& s_new, ImpulseSplitDirection& d) {
//   d.dlmd() = s_new.lmd - s.lmd;
//   d.dgmm() = s_new.gmm - s.gmm;
//   d.dmu() = s_new.mu_stack() - s.mu_stack();
//   d.df() = s_new.f_stack() - s.f_stack();
//   robot.subtractConfiguration(s_new.q, s.q, d.dq());
//   d.dv() = s_new.v - s.v;
//   impulse_dynamics_.computeCondensedDirection(kkt_matrix_, kkt_residual_, d);
//   // constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
// }


inline void SplitImpulseParNMPC::computeCondensedPrimalDirection(
    Robot& robot, const ImpulseSplitSolution& s, ImpulseSplitDirection& d) {
  impulse_dynamics_.computeCondensedPrimalDirection(robot, d);
  constraints_->computeSlackAndDualDirection(robot, constraints_data_, s, d);
}


inline void SplitImpulseParNMPC::computeCondensedDualDirection(
    const Robot& robot, const ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, ImpulseSplitDirection& d) {
  impulse_dynamics_.computeCondensedDualDirection(robot, kkt_matrix,
                                                  kkt_residual, d.dgmm(), d);
}

 
inline double SplitImpulseParNMPC::maxPrimalStepSize() {
  return constraints_->maxSlackStepSize(constraints_data_);
}


inline double SplitImpulseParNMPC::maxDualStepSize() {
  return constraints_->maxDualStepSize(constraints_data_);
}


inline void SplitImpulseParNMPC::updateDual(const double dual_step_size) {
  assert(step_size > 0);
  assert(step_size <= 1);
  constraints_->updateDual(constraints_data_, dual_step_size);
}


inline void SplitImpulseParNMPC::updatePrimal(
    const Robot& robot, const double primal_step_size, 
    const ImpulseSplitDirection& d, ImpulseSplitSolution& s, 
    const bool is_state_constraint_valid) {
  assert(primal_step_size > 0);
  assert(primal_step_size <= 1);
  s.integrate(robot, primal_step_size, d, is_state_constraint_valid);
  constraints_->updateSlack(constraints_data_, primal_step_size);
}


// void SplitImpulseParNMPC::updatePrimal(Robot& robot, const double step_size, 
//                                        const ImpulseSplitDirection& d, 
//                                        ImpulseSplitSolution& s) {
//   assert(step_size > 0);
//   assert(step_size <= 1);
//   s.lmd.noalias() += step_size * d.dlmd();
//   s.gmm.noalias() += step_size * d.dgmm();
//   s.mu_stack().noalias() += step_size * d.dmu();
//   s.set_mu_contact();
//   s.dv.noalias() += step_size * d.ddv;
//   s.f_stack().noalias() += step_size * d.df();
//   s.set_f();
//   robot.integrateConfiguration(d.dq(), step_size, s.q);
//   s.v.noalias() += step_size * d.dv();
//   s.beta.noalias() += step_size * d.dbeta;
//   // constraints_->updateSlack(constraints_data_, step_size);
// }


inline void SplitImpulseParNMPC::computeKKTResidual(
    Robot& robot, const ImpulseStatus& impulse_status, const double t,  
    const Eigen::VectorXd& q_prev, const Eigen::VectorXd& v_prev, 
    const ImpulseSplitSolution& s, const SplitSolution& s_next,
    ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual,
    const bool is_state_constraint_valid) {
  assert(q_prev.size() == robot.dimq());
  assert(v_prev.size() == robot.dimv());
  kkt_matrix.setImpulseStatus(impulse_status);
  kkt_residual.setImpulseStatus(impulse_status);
  robot.updateKinematics(s.q, s.v);
  // condensing the impulse dynamics
  kkt_residual.setZero();
  kkt_matrix.setZero();
  cost_->computeStageCostDerivatives(robot, cost_data_, t, s, kkt_residual);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, dtau, s);
  constraints_->augmentDualResidual(robot, constraints_data_, s, kkt_residual);
  stateequation::LinearizeImpulseBackwardEuler(robot, q_prev, v_prev, s, s_next, 
                                               kkt_matrix, kkt_residual);
  impulse_dynamics_.linearizeImpulseDynamics(robot, impulse_status, s,
                                             kkt_matrix, kkt_residual,
                                             is_state_constraint_valid);
}


inline double SplitImpulseParNMPC::squaredNormKKTResidual(
    const ImpulseSplitKKTResidual& kkt_residual, 
    const bool is_state_constraint_valid) const {
  double error = 0;
  error += kkt_residual.lx().squaredNorm();
  error += kkt_residual.ldv.squaredNorm();
  error += kkt_residual.lf().squaredNorm();
  error += stateequation::SquaredNormStateEuqationResidual(kkt_residual);
  error += impulse_dynamics_.squaredNormImpulseDynamicsResidual(
      kkt_residual, is_state_constraint_valid);
  error += constraints_->squaredNormPrimalAndDualResidual(constraints_data_);
  return error;
}


inline double SplitImpulseParNMPC::stageCost(Robot& robot, const double t, 
                                             const ImpulseSplitSolution& s, 
                                             const double primal_step_size) {
  assert(primal_step_size >= 0);
  assert(primal_step_size <= 1);
  robot.updateKinematics(s.q, s.v+s.dv);
  double cost = 0;
  cost += cost_->l(robot, cost_data_, t, s);
  if (primal_step_size > 0) {
    cost += constraints_->costSlackBarrier(constraints_data_, primal_step_size);
  }
  else {
    cost += constraints_->costSlackBarrier(constraints_data_);
  }
  return cost;
}


inline double SplitImpulseParNMPC::constraintViolation(
    Robot& robot, const ImpulseStatus& impulse_status, const double t, 
    const ImpulseSplitSolution& s, const Eigen::VectorXd& q_next, 
    const Eigen::VectorXd& v_next, ImpulseSplitKKTResidual& kkt_residual,
    const bool is_state_constraint_valid) {
  kkt_residual.setImpulseStatus(impulse_status);
  kkt_residual.setZero();
  robot.updateKinematics(s.q, s.v+s.dv);
  constraints_->computePrimalAndDualResidual(robot, constraints_data_, s);
  stateequation::ComputeImpulseBackwardEulerResidual(robot, s, q_next, v_next, 
                                                     kkt_residual);
  impulse_dynamics_.computeImpulseDynamicsResidual(robot, impulse_status, s, 
                                                   kkt_residual,
                                                   is_state_constraint_valid);
  double violation = 0;
  violation += constraints_->l1NormPrimalResidual(constraints_data_);
  violation += stateequation::L1NormStateEuqationResidual(kkt_residual);
  violation += impulse_dynamics_.l1NormImpulseDynamicsResidual(
      kkt_residual, is_state_constraint_valid);
  return violation;
}

} // namespace idocp