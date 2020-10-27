#ifndef IDOCP_CONDENSED_CONTACT_KKT_HXX_
#define IDOCP_CONDENSED_CONTACT_KKT_HXX_

#include "idocp/ocp/condensed_contact_kkt.hpp"

#include <assert.h>

namespace idocp {


inline void CondensedContactKKT::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  has_active_contacts_ = contact_status.hasActiveContacts();
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::dIDCdqv_() {
  return dIDCdqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::dIDdq_() {
  return dIDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::dIDdv_() {
  return dIDCdqv_full_.block(0, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::dCdq_() {
  return dIDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::dCdv_() {
  return dIDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::dCda_() {
  return dCda_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::MJtJinv_() {
  return MJtJinv_full_.topLeftCorner(dimv_+dimf_, dimv_+dimf_);
}

inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::MJtJinv_dIDCdqv_() {
  return MJtJinv_dIDCdqv_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::Qafqv_condensed_() {
  return Qafqv_condensed_full_.topLeftCorner(dimv_+dimf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> CondensedContactKKT::Qafu_condensed_() {
  return Qafu_condensed_full_.topLeftCorner(dimv_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> CondensedContactKKT::IDC_() {
  return IDC_full_.head(dimv_+dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
CondensedContactKKT::IDC_() const {
  return IDC_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> CondensedContactKKT::ID_() {
  return IDC_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
CondensedContactKKT::ID_() const {
  return IDC_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> CondensedContactKKT::C_() {
  return IDC_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
CondensedContactKKT::C_() const {
  return IDC_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
CondensedContactKKT::MJtJinv_IDC_() {
  return MJtJinv_IDC_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
CondensedContactKKT::laf_condensed_() {
  return laf_condensed_full_.head(dimv_+dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
CondensedContactKKT::la_condensed_() {
  return laf_condensed_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
CondensedContactKKT::lf_condensed_() {
  return laf_condensed_full_.segment(dimv_, dimf_);
}

} // namespace idocp 

#endif // IDOCP_CONDENSED_CONTACT_KKT_HXX_ 