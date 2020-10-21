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
                                 const SplitDirection& d_next,
                                 ImpulseSplitDirection& d);

  void computeImpulseDynamicsResidual(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const ImpulseSplitSolution& s, 
                                      ImpulseKKTResidual& kkt_residual);

  static double l1NormImpulseDynamicsResidual(
      const ImpulseKKTResidual& kkt_residual);

  static double squaredNormImpulseDynamicsResidual(
      const ImpulseKKTResidual& kkt_residual);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SchurComplement schur_complement_;
  Eigen::MatrixXd dImD_dq_, dImD_ddv_, dImD_df_full_, 
                  MJTJinv_full_, MJTJinvImDCqv_full_, 
                  Qdvq_condensed_, Qdvv_condensed_, 
                  Qfq_condensed_full_, Qfv_condensed_full_;
  Eigen::VectorXd MJTJinvImDC_full_, ldv_condensed_, lf_condensed_full_;
  int dimv_, dimf_;

  void linearizeInverseImpulseDynamics(Robot& robot, 
                                       const ContactStatus& contact_status,
                                       const ImpulseSplitSolution& s, 
                                       ImpulseKKTResidual& kkt_residual);

  static void linearizeContactConstraint(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         ImpulseKKTMatrix& kkt_matrix, 
                                         ImpulseKKTResidual& kkt_residual);

  static void setContactForces(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const ImpulseSplitSolution& s);

  static void computeInverseImpulseDynamicsResidual(
      Robot& robot, const ImpulseSplitSolution& s, 
      ImpulseKKTResidual& kkt_residual);

  static void computeContactConstraintResidual(
      const Robot& robot, const ContactStatus& contact_status, 
      ImpulseKKTResidual& kkt_residual);


  void setContactStatus(const ContactStatus& contact_status);

  Eigen::Block<Eigen::MatrixXd> dImD_df_();

  Eigen::Block<Eigen::MatrixXd> MJTJinv_();

  Eigen::Block<Eigen::MatrixXd> MJTJinvImDCqv_();

  Eigen::Block<Eigen::MatrixXd> Qfq_condensed_();

  Eigen::Block<Eigen::MatrixXd> Qfv_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> MJTJinvImDC_();

  Eigen::VectorBlock<Eigen::VectorXd> lf_condensed_();

  const Eigen::Block<const Eigen::MatrixXd> dImD_df_() const;

  const Eigen::Block<const Eigen::MatrixXd> MJTJinv_() const;

  const Eigen::Block<const Eigen::MatrixXd> MJTJinvImDCqv_() const;

  const Eigen::Block<const Eigen::MatrixXd> Qfq_condensed_() const;

  const Eigen::Block<const Eigen::MatrixXd> Qfv_condensed_() const;

  const Eigen::VectorBlock<const Eigen::VectorXd> MJTJinvImDC_() const;

  const Eigen::VectorBlock<const Eigen::VectorXd> lf_condensed_() const;

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_forward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_ 