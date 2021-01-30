#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_dynamics_forward_euler_data.hpp"


namespace idocp {

class ImpulseDynamicsForwardEuler {
public:
  ImpulseDynamicsForwardEuler(const Robot& robot, const double penalty=1.0e04);

  ImpulseDynamicsForwardEuler();

  ~ImpulseDynamicsForwardEuler();

  ImpulseDynamicsForwardEuler(const ImpulseDynamicsForwardEuler&) = default;

  ImpulseDynamicsForwardEuler& operator=(const ImpulseDynamicsForwardEuler&) 
      = default;
 
  ImpulseDynamicsForwardEuler(ImpulseDynamicsForwardEuler&&) noexcept = default;

  ImpulseDynamicsForwardEuler& operator=(ImpulseDynamicsForwardEuler&&) noexcept 
      = default;

  void linearizeImpulseDynamics(Robot& robot, 
                                const ImpulseStatus& impulse_status, 
                                const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual, 
                                const bool is_state_constraint_valid);

  static void linearizeInverseImpulseDynamics(
      Robot& robot, const ImpulseStatus& impulse_status, 
      const ImpulseSplitSolution& s, ImpulseDynamicsForwardEulerData& data);

  static void linearizeImpulseVelocityConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseDynamicsForwardEulerData& data);

  static void linearizeImpulsePositionConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, 
      ImpulseSplitKKTMatrix& kkt_matrix, ImpulseSplitKKTResidual& kkt_residual);

  void condenseImpulseDynamics(Robot& robot, 
                               const ImpulseStatus& impulse_status,
                               ImpulseSplitKKTMatrix& kkt_matrix, 
                               ImpulseSplitKKTResidual& kkt_residual);

  static void condensing(const Robot& robot, 
                         const ImpulseStatus& impulse_status,
                         ImpulseDynamicsForwardEulerData& data, 
                         ImpulseSplitKKTMatrix& kkt_matrix, 
                         ImpulseSplitKKTResidual& kkt_residual);

  void computeCondensedPrimalDirection(const Robot& robot, 
                                       ImpulseSplitDirection& d);

  template <typename VectorType>
  void computeCondensedDualDirection(const Robot& robot, 
                                     const ImpulseSplitKKTMatrix& kkt_matrix, 
                                     const ImpulseSplitKKTResidual& kkt_residual, 
                                     const Eigen::MatrixBase<VectorType>& dgmm,
                                     ImpulseSplitDirection& d);

  static void expansionPrimal(const Robot& robot, 
                              const ImpulseDynamicsForwardEulerData& data, 
                              ImpulseSplitDirection& d);

  template <typename VectorType>
  static void expansionDual(const Robot& robot, 
                            ImpulseDynamicsForwardEulerData& data,
                            const ImpulseSplitKKTMatrix& kkt_matrix, 
                            const ImpulseSplitKKTResidual& kkt_residual,
                            const Eigen::MatrixBase<VectorType>& dgmm,
                            ImpulseSplitDirection& d);

  void computeImpulseDynamicsResidual(Robot& robot, 
                                      const ImpulseStatus& impulse_status,
                                      const ImpulseSplitSolution& s, 
                                      ImpulseSplitKKTResidual& kkt_residual,
                                      const bool is_state_constraint_valid);

  double l1NormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual, 
      const bool is_state_constraint_valid) const;

  double squaredNormImpulseDynamicsResidual(
      const ImpulseSplitKKTResidual& kkt_residual,
      const bool is_state_constraint_valid) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ImpulseDynamicsForwardEulerData data_;
  double penalty_;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_forward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_ 