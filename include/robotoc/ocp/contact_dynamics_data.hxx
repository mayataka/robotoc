#ifndef ROBOTOC_CONTACT_DYNAMICS_DATA_HXX_
#define ROBOTOC_CONTACT_DYNAMICS_DATA_HXX_

#include "robotoc/ocp/contact_dynamics_data.hpp"

namespace robotoc {

inline ContactDynamicsData::ContactDynamicsData(const Robot& robot) 
  : Qxu_passive(Eigen::MatrixXd::Zero(2*robot.dimv(), robot.dim_passive())),
    Quu_passive_topRight(Eigen::MatrixXd::Zero(robot.dim_passive(), robot.dimu())),
    lu_passive(Eigen::VectorXd::Zero(robot.dim_passive())),
    dIDda(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dCda_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    dIDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        2*robot.dimv())),
    MJtJinv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                        robot.dimv()+robot.max_dimf())), 
    MJtJinv_dIDCdqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                                2*robot.dimv())), 
    Qafqv_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                      2*robot.dimv())), 
    Qafu_full_full_(Eigen::MatrixXd::Zero(robot.dimv()+robot.max_dimf(), 
                                          robot.dimv())), 
    IDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    MJtJinv_IDC_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    laf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    haf_full_(Eigen::VectorXd::Zero(robot.dimv()+robot.max_dimf())),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    dimf_(0),
    dimvf_(robot.dimv()),
    dim_passive_(robot.dim_passive()) {
}


inline ContactDynamicsData::ContactDynamicsData() 
  : Qxu_passive(),
    Quu_passive_topRight(),
    lu_passive(),
    dIDda(),
    dIDCdqv_full_(),
    MJtJinv_full_(), 
    MJtJinv_dIDCdqv_full_(), 
    Qafqv_full_(), 
    Qafu_full_full_(), 
    IDC_full_(),
    MJtJinv_IDC_full_(),
    laf_full_(),
    haf_full_(),
    dimv_(0),
    dimu_(0),
    dimf_(0),
    dimvf_(0),
    dim_passive_(0) {
}


inline ContactDynamicsData::~ContactDynamicsData() {
}


inline void ContactDynamicsData::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  dimvf_ = dimv_ + dimf_;
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::dCda() {
  return dCda_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::dCda() const {
  return dCda_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::dIDCdqv() {
  return dIDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::dIDCdqv() const {
  return dIDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::dIDdq() {
  return dIDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::dIDdq() const {
  return dIDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::dIDdv() {
  return dIDCdqv_full_.block(0, dimv_, dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::dIDdv() const {
  return dIDCdqv_full_.block(0, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::dCdq() {
  return dIDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::dCdq() const {
  return dIDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::dCdv() {
  return dIDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::dCdv() const {
  return dIDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::MJtJinv() {
  return MJtJinv_full_.topLeftCorner(dimvf_, dimvf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::MJtJinv() const {
  return MJtJinv_full_.topLeftCorner(dimvf_, dimvf_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::MJtJinv_dIDCdqv() {
  return MJtJinv_dIDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::MJtJinv_dIDCdqv() const {
  return MJtJinv_dIDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::Qafqv() {
  return Qafqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::Qafqv() const {
  return Qafqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::Qafu_full() {
  return Qafu_full_full_.topLeftCorner(dimvf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::Qafu_full() const {
  return Qafu_full_full_.topLeftCorner(dimvf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::Qafu_passive() {
  return Qafu_full_full_.topLeftCorner(dimvf_, dim_passive_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::Qafu_passive() const {
  return Qafu_full_full_.topLeftCorner(dimvf_, dim_passive_);
}


inline Eigen::Block<Eigen::MatrixXd> ContactDynamicsData::Qafu() {
  return Qafu_full_full_.block(0, dim_passive_, dimvf_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ContactDynamicsData::Qafu() const {
  return Qafu_full_full_.block(0, dim_passive_, dimvf_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::IDC() {
  return IDC_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::IDC() const {
  return IDC_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::ID_full() {
  return IDC_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::ID_full() const {
  return IDC_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::ID_passive() {
  return IDC_full_.head(dim_passive_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::ID_passive() const {
  return IDC_full_.head(dim_passive_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::ID() {
  return IDC_full_.segment(dim_passive_, dimu_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::ID() const {
  return IDC_full_.segment(dim_passive_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::C() {
  return IDC_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::C() const {
  return IDC_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::MJtJinv_IDC() {
  return MJtJinv_IDC_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::MJtJinv_IDC() const {
  return MJtJinv_IDC_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::laf() {
  return laf_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::laf() const {
  return laf_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::la() {
  return laf_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::la() const {
  return laf_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::lf() {
  return laf_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::lf() const {
  return laf_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::haf() {
  return haf_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::haf() const {
  return haf_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::ha() {
  return haf_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::ha() const {
  return haf_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> ContactDynamicsData::hf() {
  return haf_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ContactDynamicsData::hf() const {
  return haf_full_.segment(dimv_, dimf_);
}

} // namespace robotoc 

#endif // ROBOTOC_CONTACT_DYNAMICS_DATA_HXX_ 