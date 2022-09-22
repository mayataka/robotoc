#ifndef ROBOTOC_IMPULSE_DYNAMICS_DATA_HXX_
#define ROBOTOC_IMPULSE_DYNAMICS_DATA_HXX_

#include "robotoc/dynamics/impulse_dynamics_data.hpp"

namespace robotoc {

inline void ImpulseDynamicsData::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimf();
  dimvf_ = dimv_ + dimf_;
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsData::dImDCdqv() {
  return dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::dImDCdqv() const {
  return dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsData::dImDCdq() {
  return dImDCdqv_full_.topLeftCorner(dimvf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::dImDCdq() const {
  return dImDCdqv_full_.topLeftCorner(dimvf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsData::dImDdq() {
  return dImDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::dImDdq() const {
  return dImDCdqv_full_.topLeftCorner(dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsData::dCdq() {
  return dImDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::dCdq() const {
  return dImDCdqv_full_.block(dimv_, 0, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsData::dCdv() {
  return dImDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::dCdv() const {
  return dImDCdqv_full_.block(dimv_, dimv_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsData::MJtJinv() {
  return MJtJinv_full_.topLeftCorner(dimvf_, dimvf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::MJtJinv() const {
  return MJtJinv_full_.topLeftCorner(dimvf_, dimvf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseDynamicsData::MJtJinv_dImDCdqv() {
  return MJtJinv_dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::MJtJinv_dImDCdqv() const {
  return MJtJinv_dImDCdqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseDynamicsData::Qdvfqv() {
  return Qdvfqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseDynamicsData::Qdvfqv() const {
  return Qdvfqv_full_.topLeftCorner(dimvf_, 2*dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsData::ImDC() {
  return ImDC_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsData::ImDC() const {
  return ImDC_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsData::ImD() {
  return ImDC_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsData::ImD() const {
  return ImDC_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsData::C() {
  return ImDC_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsData::C() const {
  return ImDC_full_.segment(dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsData::MJtJinv_ImDC() {
  return MJtJinv_ImDC_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsData::MJtJinv_ImDC() const {
  return MJtJinv_ImDC_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsData::ldvf() {
  return ldvf_full_.head(dimvf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsData::ldvf() const {
  return ldvf_full_.head(dimvf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsData::ldv() {
  return ldvf_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsData::ldv() const {
  return ldvf_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseDynamicsData::lf() {
  return ldvf_full_.segment(dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseDynamicsData::lf() const {
  return ldvf_full_.segment(dimv_, dimf_);
}

} // namespace robotoc 

#endif // ROBOTOC_IMPULSE_DYNAMICS_DATA_HXX_ 