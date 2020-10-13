#ifndef IDOCP_IMPULSE_DYNAMICS_HPP_
#define IDOCP_IMPULSE_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"


namespace idocp {

class ImpulseDynamics {
public:
  ImpulseDynamics(const Robot& robot);

  ImpulseDynamics();

  ~ImpulseDynamics();

  ImpulseDynamics(const ImpulseDynamics&) = default;

  ImpulseDynamics& operator=(const ImpulseDynamics&) = default;
 
  ImpulseDynamics(ImpulseDynamics&&) noexcept = default;

  ImpulseDynamics& operator=(ImpulseDynamics&&) noexcept = default;

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
  Eigen::MatrixXd Qdvdv_, dImD_dq_, dImD_ddv_, dImD_df_full_, 
                  dC_ddv_;
  int dimf_;

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

  const Eigen::Block<const Eigen::MatrixXd, Eigen::Dynamic, Eigen::Dynamic, true> 
  dImD_df_() const;

};

} // namespace idocp 

#include "idocp/impulse/impulse_dynamics.hxx"

#endif // IDOCP_IMPULSE_DYNAMICS_HPP_ 