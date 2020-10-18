#ifndef IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_
#define IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/ocp/schur_complement.hpp"


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
                                const ContactStatus& contact_status, 
                                const ImpulseSplitSolution& s, 
                                ImpulseKKTMatrix& kkt_matrix, 
                                ImpulseKKTResidual& kkt_residual);

  void condenseImpulseDynamics(Robot& robot, 
                               const ContactStatus& contact_status,
                               const ImpulseSplitSolution& s, 
                               ImpulseKKTMatrix& kkt_matrix, 
                               ImpulseKKTResidual& kkt_residual);

  void computeCondensedDirection(const ImpulseKKTMatrix& kkt_matrix, 
                                 const ImpulseKKTResidual& kkt_residual, 
                                 const ImpulseSplitDirection& d_next,
                                 ImpulseSplitDirection& d);

  void computeImpulseDynamicsResidual(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const ImpulseSplitSolution& s, 
                                      ImpulseKKTResidual& kkt_residual);

  double l1NormImpulseDynamicsResidual(
      const ImpulseKKTResidual& kkt_residual) const;

  double squaredNormImpulseDynamicsResidual(
      const ImpulseKKTResidual& kkt_residual) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SchurComplement schur_complement_;
  Eigen::MatrixXd dImD_dq_, dImD_ddv_, dImD_df_full_, 
                  MJTJinv_full_, MJTJinvImDCqv_full_,
                  Qdvq_, Qdvv_, Qfq_full_, Qfv_full_;
  Eigen::VectorXd ldvf_full_, ImDC_full_, ldvf_condensed_full_;
  int dimv_, dimf_;

  void linearizeInverseImpulseDynamics(Robot& robot, 
                                       const ContactStatus& contact_status,
                                       const ImpulseSplitSolution& s, 
                                       ImpulseKKTResidual& kkt_residual);

  static void linearizeContactVelocityConstraint(
      Robot& robot, const ContactStatus& contact_status, 
      ImpulseKKTMatrix& kkt_matrix, ImpulseKKTResidual& kkt_residual);

  static void setContactForces(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const ImpulseSplitSolution& s);

  static void computeInverseImpulseDynamicsResidual(
      Robot& robot, const ImpulseSplitSolution& s, 
      ImpulseKKTResidual& kkt_residual);

  static void computeContactVelocityConstraintResidual(
      const Robot& robot, const ContactStatus& contact_status, 
      ImpulseKKTResidual& kkt_residual);


  void setContactStatus(const ContactStatus& contact_status);

  Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  dImD_df_();

  Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  MJTJinv_();

  Eigen::Block<Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  MJTJinvImDCqv_();

  Eigen::VectorBlock<Eigen::VectorXd> ldvf_();

  Eigen::VectorBlock<Eigen::VectorXd> ImDC_();

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_forward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_ 