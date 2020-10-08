#ifndef IDOCP_IMPULSE_DYNAMICS_HPP_
#define IDOCP_IMPULSE_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {
namespace impulsedynamics {
  
void LinearizeImpulseDynamics(Robot& robot, 
                              const ContactStatus& contact_status, 
                              const SplitSolution& s, KKTMatrix& kkt_matrix, 
                              KKTResidual& kkt_residual);

void ComputeImpulseDynamicsResidual(Robot& robot, 
                                    const ContactStatus& contact_status,
                                    const SplitSolution& s, 
                                    KKTResidual& kkt_residual);

double L1NormImpulseDynamicsResidual(const KKTResidual& kkt_residual);

double SquaredNormImpulseDynamicsResidual(const KKTResidual& kkt_residual);

void LinearizeInverseImpulseDynamics(Robot& robot, 
                                     const ContactStatus& contact_status,
                                     const SplitSolution& s, 
                                     KKTResidual& kkt_residual);

void LinearizeContactConstraint(Robot& robot,   
                                const ContactStatus& contact_status, 
                                KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual);

void SetContactForces(Robot& robot, const ContactStatus& contact_status, 
                      const SplitSolution& s);

void ComputeInverseDynamicsResidual(Robot& robot, const SplitSolution& s, 
                                    KKTResidual& kkt_residual);

void ComputeFloatingBaseConstraintResidual(const Robot& robot, 
                                           const SplitSolution& s, 
                                           KKTResidual& kkt_residual);

void ComputeContactConstraintResidual(const Robot& robot, 
                                      const ContactStatus& contact_status, 
                                      KKTResidual& kkt_residual);

} // namespace impulsedynamics
} // namespace idocp 

#include "idocp/ocp/impulse_dynamics.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_HPP_ 