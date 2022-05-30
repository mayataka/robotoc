#ifndef ROBOTOC_SPLIT_OCP_HXX_
#define ROBOTOC_SPLIT_OCP_HXX_

#include "robotoc/ocp/split_ocp.hpp"

#include <cassert>

namespace robotoc {

template <typename SplitSolutionType>
inline void SplitOCP::computeKKTResidual_impl(Robot& robot, 
                                              const ContactStatus& contact_status, 
                                              const GridInfo& grid_info, 
                                              const Eigen::VectorXd& q_prev, 
                                              const SplitSolution& s,
                                              const SplitSolutionType& s_next, 
                                              SplitKKTMatrix& kkt_matrix,
                                              SplitKKTResidual& kkt_residual) {
  robot.updateKinematics(s.q, s.v, s.a);
  kkt_matrix.setContactStatus(contact_status);
  kkt_residual.setContactStatus(contact_status);
  kkt_residual.setZero();
  stage_cost_ = cost_->linearizeStageCost(robot, contact_status, cost_data_,  
                                          grid_info, s, kkt_residual);
  kkt_residual.h = (1.0/grid_info.dt) * stage_cost_;
  constraints_->linearizeConstraints(robot, contact_status, constraints_data_, 
                                     s, kkt_residual);
  barrier_cost_ = constraints_data_.logBarrier();
  state_equation_.linearizeStateEquation(robot, grid_info.dt, q_prev, s, s_next, 
                                         kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, s, 
                                             kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
}


template <typename SplitSolutionType>
inline void SplitOCP::computeKKTSystem_impl(Robot& robot, 
                                            const ContactStatus& contact_status,  
                                            const GridInfo& grid_info, 
                                            const Eigen::VectorXd& q_prev, 
                                            const SplitSolution& s, 
                                            const SplitSolutionType& s_next,
                                            SplitKKTMatrix& kkt_matrix, 
                                            SplitKKTResidual& kkt_residual) {
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
  barrier_cost_ = constraints_data_.logBarrier();
  state_equation_.linearizeStateEquation(robot, grid_info.dt, q_prev, s, s_next, 
                                         kkt_matrix, kkt_residual);
  contact_dynamics_.linearizeContactDynamics(robot, contact_status, s,
                                             kkt_residual);
  kkt_residual.kkt_error = KKTError(kkt_residual);
  constraints_->condenseSlackAndDual(contact_status, constraints_data_, 
                                     kkt_matrix, kkt_residual);
  contact_dynamics_.condenseContactDynamics(robot, contact_status, grid_info.dt, 
                                            kkt_matrix, kkt_residual);
  state_equation_.correctLinearizedStateEquation(robot, grid_info.dt, s, s_next, 
                                                 kkt_matrix, kkt_residual);
}


template <typename SplitDirectionType>
inline void SplitOCP::expandDual_impl(const GridInfo& grid_info, 
                                      const SplitDirectionType& d_next, 
                                      SplitDirection& d, const double dts) {
  assert(grid_info.dt > 0);
  contact_dynamics_.expandDual(grid_info.dt, dts, d_next, d);
  state_equation_.correctCostateDirection(d);
}

} // namespace robotoc

#endif // ROBOTOC_SPLIT_OCP_HXX_