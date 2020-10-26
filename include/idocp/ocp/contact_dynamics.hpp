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

  void linearizeContactDynamics(Robot& robot, 
                                const ContactStatus& contact_status, 
                                const double dtau, const SplitSolution& s, 
                                KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual);

  void condenseContactDynamics(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const double dtau, const SplitSolution& s, 
                               KKTMatrix& kkt_matrix, 
                               KKTResidual& kkt_residual);

  template <typename VectorType>
  void computeCondensedDirection(const double dtau, const SplitSolution& s, 
                                 const KKTMatrix& kkt_matrix, 
                                 const KKTResidual& kkt_residual, 
                                 const Eigen::MatrixBase<VectorType>& dgmm,
                                 SplitDirection& d);

  void computeContactDynamicsResidual(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const double dtau, const SplitSolution& s, 
                                      KKTResidual& kkt_residual);

  double l1NormContactDynamicsResidual(const double dtau) const;

  double squaredNormContactDynamicsResidual(const double dtau) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd dIDda_, dCda_full_, dIDCdqv_full_, MJtJinv_full_, 
                  MJtJinv_dIDCdqv_full_, Qafqv_condensed_full_, 
                  Qafu_condensed_full_;
  Eigen::VectorXd IDC_full_, MJtJinv_IDC_full_, laf_condensed_full_;
  int dimv_, dimu_, dimf_, dim_passive_;
  bool has_floating_base_, has_active_contacts_;
  static constexpr int kDimFloatingBase = 6;

  void linearizeInverseDynamics(Robot& robot, 
                                const ContactStatus& contact_status,
                                const SplitSolution& s, 
                                KKTResidual& kkt_residual);

  void linearizeContactConstraint(Robot& robot, 
                                  const ContactStatus& contact_status, 
                                  const double dtau, KKTMatrix& kkt_matrix, 
                                  KKTResidual& kkt_residual);

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

  Eigen::VectorBlock<Eigen::VectorXd> IDC_();

  const Eigen::VectorBlock<const Eigen::VectorXd> IDC_() const;

  Eigen::VectorBlock<Eigen::VectorXd> ID_();

  Eigen::VectorBlock<Eigen::VectorXd> C_();

  Eigen::VectorBlock<Eigen::VectorXd> MJtJinv_IDC_();

  Eigen::VectorBlock<Eigen::VectorXd> laf_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> la_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> lf_condensed_();

};

} // namespace idocp 

#include "idocp/ocp/contact_dynamics.hxx"

#endif // IDOCP_CONTACT_DYNAMICS_HPP_ 