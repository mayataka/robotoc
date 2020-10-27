#ifndef IDOCP_CONDENSED_CONTACT_KKT_HPP_
#define IDOCP_CONDENSED_CONTACT_KKT_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


namespace idocp {

class CondensedContactKKT {
public:
  CondensedContactKKT(const Robot& robot);

  CondensedContactKKT();

  ~CondensedContactKKT();

  CondensedContactKKT(const CondensedContactKKT&) = default;

  CondensedContactKKT& operator=(const CondensedContactKKT&) 
      = default;
 
  CondensedContactKKT(CondensedContactKKT&&) noexcept = default;

  CondensedContactKKT& operator=(CondensedContactKKT&&) noexcept 
      = default;

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

  const Eigen::VectorBlock<const Eigen::VectorXd> ID_() const;

  Eigen::VectorBlock<Eigen::VectorXd> C_();

  const Eigen::VectorBlock<const Eigen::VectorXd> C_() const;

  Eigen::VectorBlock<Eigen::VectorXd> MJtJinv_IDC_();

  Eigen::VectorBlock<Eigen::VectorXd> laf_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> la_condensed_();

  Eigen::VectorBlock<Eigen::VectorXd> lf_condensed_();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd dIDda_, dCda_full_, dIDCdqv_full_, MJtJinv_full_, 
                  MJtJinv_dIDCdqv_full_, Qafqv_condensed_full_, 
                  Qafu_condensed_full_;
  Eigen::VectorXd u_passive_, IDC_full_, MJtJinv_IDC_full_, laf_condensed_full_;
  int dimv_, dimu_, dimf_, dim_passive_;
  bool has_floating_base_, has_active_contacts_;
  static constexpr int kDimFloatingBase = 6;

};

} // namespace idocp 

#include "idocp/ocp/condensed_contact_kkt.hxx"

#endif // IDOCP_CONDENSED_CONTACT_KKT_HPP_ 