#ifndef IDOCP_CONTACT_DYNAMICS_HPP_
#define IDOCP_CONTACT_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class ContactDynamics {
public:
  ContactDynamics(const Robot& robot);

  ContactDynamics();

  ~ContactDynamics();

  ContactDynamics(const ContactDynamics&) = default;

  ContactDynamics& operator=(const ContactDynamics&) 
      = default;
 
  ContactDynamics(ContactDynamics&&) noexcept = default;

  ContactDynamics& operator=(ContactDynamics&&) noexcept 
      = default;

  void linearizeRobotDynamics(Robot& robot, const ContactStatus& contact_status, 
                              const double dtau, const SplitSolution& s, 
                              KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  void condenseRobotDynamics(Robot& robot, const ContactStatus& contact_status,
                             const double dtau, const SplitSolution& s, 
                             KKTMatrix& kkt_matrix, KKTResidual& kkt_residual);

  void computeCondensedDirection(const double dtau, const KKTMatrix& kkt_matrix, 
                                 const KKTResidual& kkt_residual, 
                                 const SplitDirection& d_next,
                                 SplitDirection& d);

  void computeRobotDynamicsResidual(Robot& robot, 
                                    const ContactStatus& contact_status,
                                    const double dtau, const SplitSolution& s, 
                                    KKTResidual& kkt_residual);

  static double l1NormRobotDynamicsResidual(const double dtau, 
                                            const KKTResidual& kkt_residual);

  static double squaredNormRobotDynamicsResidual(
      const double dtau, const KKTResidual& kkt_residual);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd dIDC_dqv_full_, dIDC_daf_full_, MJtJinv_full_, 
                  MJtJinv_IDCqv_full_, Qafqv_condensed_full_, 
                  Qafu_condensed_full_;
  Eigen::VectorXd MJtJinv_IDC_full_, laf_condensed_full_;
  int dimv_, dimf_;

  void linearizeInverseDynamics(Robot& robot, 
                                const ContactStatus& contact_status,
                                const double dtau, const SplitSolution& s, 
                                KKTResidual& kkt_residual);

  static void linearizeContactConstraint(Robot& robot, 
                                         const ContactStatus& contact_status, 
                                         const double dtau,
                                         KKTMatrix& kkt_matrix, 
                                         KKTResidual& kkt_residual);

  static void setContactForces(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const SplitSolution& s);

  static void computeInverseDynamicsResidual(
      Robot& robot, const double dtau, const SplitSolution& s, 
      KKTResidual& kkt_residual);

  void computeFloatingBaseConstraintResidual(const Robot& robot, 
                                             const double dtau,
                                             const SplitSolution& s, 
                                             KKTResidual& kkt_residual);

  static void computeContactConstraintResidual(
      const Robot& robot, const ContactStatus& contact_status, 
      const double dtau, KKTResidual& kkt_residual);


  void setContactStatus(const ContactStatus& contact_status);

  Eigen::Block<Eigen::MatrixXd> dIDCdqv_();

  Eigen::Block<Eigen::MatrixXd> dIDdq_();

  Eigen::Block<Eigen::MatrixXd> dIDdv_();

  Eigen::Block<Eigen::MatrixXd> dCdq_();

  Eigen::Block<Eigen::MatrixXd> dCdv_();

  Eigen::Block<Eigen::MatrixXd> dCda_();

  Eigen::Block<Eigen::MatrixXd> MJtJinv_();

  Eigen::Block<Eigen::MatrixXd> MJtJinv_dIDCdqv_();

  Eigen::Block<Eigen::MatrixXd> Qafqv_condensed_();

  Eigen::Block<Eigen::MatrixXd> Qafu_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> MJtJinv_IDC_();

  Eigen::VectorBlock<Eigen::VectorXd> laf_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> la_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> lf_condensed_();

};

} // namespace idocp 

#include "idocp/impulse/contact_dynamics.hxx"

#endif // IDOCP_CONTACT_DYNAMICS_HPP_ 