#ifndef ROBOTOC_SWITCHING_CONSTRAINT_JACOBIAN_HXX_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_JACOBIAN_HXX_

#include "robotoc/core/switching_constraint_jacobian.hpp"

#include <cassert>


namespace robotoc {

inline void SwitchingConstraintJacobian::setDimension(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimi();
}


inline void SwitchingConstraintJacobian::setDimension() {
  dimi_ = 0;
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintJacobian::Pq() {
  return Pq_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintJacobian::Pq() const {
  return Pq_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintJacobian::Phix() {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintJacobian::Phix() const {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintJacobian::Phiq() {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintJacobian::Phiq() const {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintJacobian::Phiv() {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintJacobian::Phiv() const {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintJacobian::Phia() {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintJacobian::Phia() const {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintJacobian::Phiu() {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintJacobian::Phiu() const {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SwitchingConstraintJacobian::Phit() {
  return Phit_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SwitchingConstraintJacobian::Phit() const {
  return Phit_full_.head(dimi_);
}


inline void SwitchingConstraintJacobian::setZero() {
  Pq().setZero();
  Phix().setZero();
  Phia().setZero();
  Phiu().setZero();
  Phit().setZero();
}


inline int SwitchingConstraintJacobian::dimi() const {
  return dimi_;
}

} // namespace robotoc 

#endif // ROBOTOC_SWITCHING_CONSTRAINT_JACOBIAN_HXX_ 