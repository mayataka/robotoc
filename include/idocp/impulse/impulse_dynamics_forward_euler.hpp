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
  Eigen::MatrixXd dImD_dq_, dImD_ddv_, MJtJinv_full_, MJtJinv_dImDCdqv_full_, 
                  Qdvfqv_condensed_full_;
  Eigen::VectorXd MJtJinv_ImDC_full_, ldvf_condensed_full_;
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

  Eigen::Block<Eigen::MatrixXd> MJtJinv_();

  Eigen::Block<Eigen::MatrixXd> MJtJinv_dImDCdqv_();

  Eigen::Block<Eigen::MatrixXd> Qdvfqv_condensed_();

  Eigen::Block<Eigen::MatrixXd> Qdvq_condensed_();

  Eigen::Block<Eigen::MatrixXd> Qdvv_condensed_();

  Eigen::Block<Eigen::MatrixXd> Qfq_condensed_();

  Eigen::Block<Eigen::MatrixXd> Qfv_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> MJtJinv_ImDC_();

  Eigen::VectorBlock<Eigen::VectorXd> ldvf_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> ldv_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> lf_condensed_();

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics_forward_euler.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_FORWARD_EULER_HPP_ 