#ifndef IDOCP_SPLIT_SWITCHING_CONSTRAINT_JACOBIAN_HXX_ 
#define IDOCP_SPLIT_SWITCHING_CONSTRAINT_JACOBIAN_HXX_

#include "idocp/ocp/split_switching_constraint_jacobian.hpp"

#include <cassert>


namespace idocp {

inline SplitSwitchingConstraintJacobian::SplitSwitchingConstraintJacobian(
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


inline SplitSwitchingConstraintJacobian::SplitSwitchingConstraintJacobian()
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


inline SplitSwitchingConstraintJacobian::~SplitSwitchingConstraintJacobian() {
}


inline void SplitSwitchingConstraintJacobian::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimf();
}


inline void SplitSwitchingConstraintJacobian::setImpulseStatus() {
  dimi_ = 0;
}


inline Eigen::Block<Eigen::MatrixXd> SplitSwitchingConstraintJacobian::Pq() {
  return Pq_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitSwitchingConstraintJacobian::Pq() const {
  return Pq_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitSwitchingConstraintJacobian::Phix() {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitSwitchingConstraintJacobian::Phix() const {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitSwitchingConstraintJacobian::Phiq() {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitSwitchingConstraintJacobian::Phiq() const {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitSwitchingConstraintJacobian::Phiv() {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitSwitchingConstraintJacobian::Phiv() const {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitSwitchingConstraintJacobian::Phia() const {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitSwitchingConstraintJacobian::Phia() {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitSwitchingConstraintJacobian::Phiu() const {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitSwitchingConstraintJacobian::Phiu() {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitSwitchingConstraintJacobian::Phit() {
  return Phit_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSwitchingConstraintJacobian::Phit() const {
  return Phit_full_.head(dimi_);
}


inline void SplitSwitchingConstraintJacobian::setZero() {
  Pq().setZero();
  Phix().setZero();
  Phia().setZero();
  Phiu().setZero();
  Phit().setZero();
}


inline int SplitSwitchingConstraintJacobian::dimi() const {
  return dimi_;
}


inline bool SplitSwitchingConstraintJacobian::isApprox(
    const SplitSwitchingConstraintJacobian& other) const {
  assert(dimi() == other.dimi());
  if (!Pq().isApprox(other.Pq())) return false;
  if (!Phix().isApprox(other.Phix())) return false;
  if (!Phia().isApprox(other.Phia())) return false;
  if (!Phiu().isApprox(other.Phiu())) return false;
  if (!Phit().isApprox(other.Phit())) return false;
  return true;
}


inline bool SplitSwitchingConstraintJacobian::hasNaN() const {
  if (Pq().hasNaN()) return true;
  if (Phix().hasNaN()) return true;
  if (Phia().hasNaN()) return true;
  if (Phiu().hasNaN()) return true;
  if (Phit().hasNaN()) return true;
  return false;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_SWITCHING_CONSTRAINT_JACOBIAN_HXX_ 