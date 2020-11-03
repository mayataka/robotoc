#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

class ImpulseDynamicsForwardEuler {
public:
  ImpulseDynamicsForwardEuler(const Robot& robot);

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
                                ImpulseKKTMatrix& kkt_matrix, 
                                ImpulseKKTResidual& kkt_residual);

  static void linearizeInverseDynamics(Robot& robot, 
                                       const ImpulseStatus& impulse_status,
                                       const ImpulseSplitSolution& s, 
                                       ImpulseDynamicsData& data);

  static void linearizeContactConstraint(Robot& robot, 
                                         const ImpulseStatus& impulse_status,
                                         ImpulseDynamicsData& data);

  void condenseImpulseDynamics(Robot& robot, 
                               const ImpulseStatus& impulse_status,
                               ImpulseKKTMatrix& kkt_matrix, 
                               ImpulseKKTResidual& kkt_residual);

  static void condensing(const Robot& robot, ImpulseDynamicsData& data, 
                         ImpulseKKTMatrix& kkt_matrix, 
                         ImpulseKKTResidual& kkt_residual);

  void computeCondensedPrimalDirection(const Robot& robot, 
                                       ImpulseSplitDirection& d);

  template <typename VectorType>
  void computeCondensedDualDirection(const Robot& robot, 
                                     const ImpulseKKTMatrix& kkt_matrix, 
                                     const ImpulseKKTResidual& kkt_residual, 
                                     const Eigen::MatrixBase<VectorType>& dgmm,
                                     ImpulseSplitDirection& d);

  static void expansionPrimal(const Robot& robot, 
                              const ImpulseDynamicsData& data, 
                              ImpulseSplitDirection& d);

  template <typename VectorType>
  static void expansionDual(const Robot& robot, ImpulseDynamicsData& data,
                            const ImpulseKKTMatrix& kkt_matrix, 
                            const ImpulseKKTResidual& kkt_residual,
                            const Eigen::MatrixBase<VectorType>& dgmm,
                            ImpulseSplitDirection& d);


  void computeImpulseDynamicsResidual(Robot& robot, 
                                      const ImpulseStatus& impulse_status,
                                      const ImpulseSplitSolution& s, 
                                      ImpulseKKTResidual& kkt_residual);

  double l1NormImpulseDynamicsResidual() const;

  double squaredNormImpulseDynamicsResidual() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  ImpulseDynamicsData data_;
  int dimv_, dimf_, dim_passive_;
  bool has_active_impulse_;

  Eigen::MatrixXd dImD_dq_, dImD_ddv_, MJtJinv_full_, MJtJinv_dImDCdqv_full_, 
                  Qdvfqv_condensed_full_;
  Eigen::VectorXd MJtJinv_ImDC_full_, ldvf_condensed_full_;
  int dimv_, dimf_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_forward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_ 