#ifndef ROBOTOC_SWITCHING_CONSTRAINT_JACOBIAN_HXX_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_JACOBIAN_HXX_

#include "robotoc/core/switching_constraint_jacobian.hpp"

#include <cassert>


namespace robotoc {

inline SwitchingConstraintJacobian::SwitchingConstraintJacobian(
    const Robot& robot)
  : Pq_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Phix_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    Phia_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Phiu_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    Phit_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimu_(robot.dimu()),
    dimi_(0) {
}


inline SwitchingConstraintJacobian::SwitchingConstraintJacobian()
  : Pq_full_(),
    Phix_full_(),
    Phia_full_(),
    Phiu_full_(),
    Phit_full_(),
    has_floating_base_(false),
    dimv_(0),
    dimx_(0),
    dimu_(0),
    dimi_(0) {
}


inline SwitchingConstraintJacobian::~SwitchingConstraintJacobian() {
}


inline void SwitchingConstraintJacobian::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimi();
}


inline void SwitchingConstraintJacobian::setImpulseStatus() {
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


inline bool SwitchingConstraintJacobian::isApprox(
    const SwitchingConstraintJacobian& other) const {
  assert(dimi() == other.dimi());
  if (!Pq().isApprox(other.Pq())) return false;
  if (!Phix().isApprox(other.Phix())) return false;
  if (!Phia().isApprox(other.Phia())) return false;
  if (!Phiu().isApprox(other.Phiu())) return false;
  if (!Phit().isApprox(other.Phit())) return false;
  return true;
}


inline bool SwitchingConstraintJacobian::hasNaN() const {
  if (Pq().hasNaN()) return true;
  if (Phix().hasNaN()) return true;
  if (Phia().hasNaN()) return true;
  if (Phiu().hasNaN()) return true;
  if (Phit().hasNaN()) return true;
  return false;
}

} // namespace robotoc 

#endif // ROBOTOC_SWITCHING_CONSTRAINT_JACOBIAN_HXX_ 