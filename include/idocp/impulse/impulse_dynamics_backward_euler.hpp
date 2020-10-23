#ifndef IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_
#define IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/ocp/schur_complement.hpp"


namespace idocp {

class ImpulseDynamicsBackwardEuler {
public:
  ImpulseDynamicsBackwardEuler(const Robot& robot);

  ImpulseDynamicsBackwardEuler();

  ~ImpulseDynamicsBackwardEuler();

  ImpulseDynamicsBackwardEuler(const ImpulseDynamicsBackwardEuler&) = default;

  ImpulseDynamicsBackwardEuler& operator=(const ImpulseDynamicsBackwardEuler&) = default;
 
  ImpulseDynamicsBackwardEuler(ImpulseDynamicsBackwardEuler&&) noexcept = default;

  ImpulseDynamicsBackwardEuler& operator=(ImpulseDynamicsBackwardEuler&&) noexcept = default;

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
                                 ImpulseSplitDirection& d) const;

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
  Eigen::MatrixXd dImD_dq_, dImD_ddv_, dImD_df_full_, Minv_, MinvImDq_, 
                  MinvImDf_full_, Qdvq_condensed_, Qdvf_condensed_full_;
  Eigen::VectorXd MinvImD_, ldv_condensed_;
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

  const Eigen::Block<const Eigen::MatrixXd> dImD_df_() const;

  Eigen::Block<Eigen::MatrixXd> MinvImDf_();

  const Eigen::Block<const Eigen::MatrixXd> MinvImDf_() const;

  Eigen::Block<Eigen::MatrixXd> Qdvf_condensed_();

  const Eigen::Block<const Eigen::MatrixXd> Qdvf_condensed_() const;

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_backward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_BACKWARD_EULER_HPP_ 