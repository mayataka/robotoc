#ifndef IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_
#define IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_backward_euler_data.hpp"


namespace idocp {

class ImpulseDynamicsBackwardEuler {
public:
  ImpulseDynamicsBackwardEuler(const Robot& robot);

  ImpulseDynamicsBackwardEuler();

  ~ImpulseDynamicsBackwardEuler();

  ImpulseDynamicsBackwardEuler(
      const ImpulseDynamicsBackwardEuler&) = default;

  ImpulseDynamicsBackwardEuler& operator=(
      const ImpulseDynamicsBackwardEuler&) = default;
 
  ImpulseDynamicsBackwardEuler(
      ImpulseDynamicsBackwardEuler&&) noexcept = default;

  ImpulseDynamicsBackwardEuler& operator=(
      ImpulseDynamicsBackwardEuler&&) noexcept = default;

  void linearizeImpulseDynamics(Robot& robot, 
                                const ImpulseStatus& impulse_status,  
                                const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual);

  static void linearizeInverseImpulseDynamics(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const ImpulseSplitSolution& s, ImpulseDynamicsBackwardEulerData& data);

  static void linearizeImpulseVelocityConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual);

  static void linearizeImpulsePositionConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual);

  void condenseImpulseDynamics(Robot& robot, 
                               const ImpulseStatus& impulse_status,  
                               ImpulseSplitKKTMatrix& kkt_matrix, 
                               ImpulseSplitKKTResidual& kkt_residual);

  static void condensing(const Robot& robot, 
                         ImpulseDynamicsBackwardEulerData& data, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual);

  void computeCondensedPrimalDirection(const Robot& robot, 
                                       const ImpulseSplitKKTMatrix& kkt_matrix, 
                                       ImpulseSplitDirection& d) const;

  void computeCondensedDualDirection(const Robot& robot, 
                                     ImpulseSplitDirection& d);

  void computeImpulseDynamicsResidual(Robot& robot, 
                                      const ImpulseStatus& impulse_status,
                                      const ImpulseSplitSolution& s, 
                                      ImpulseSplitKKTResidual& kkt_residual);

  double l1NormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

  double squaredNormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ImpulseDynamicsBackwardEulerData data_;
  int dimv_, dimf_;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_backward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_ 