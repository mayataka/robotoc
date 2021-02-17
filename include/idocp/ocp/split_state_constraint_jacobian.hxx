#ifndef IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HXX_ 
#define IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HXX_

#include "idocp/ocp/split_state_constraint_jacobian.hpp"


namespace idocp {

inline SplitStateConstraintJacobian::SplitStateConstraintJacobian(
    const Robot& robot) 
  : dintegrate_dq(),
    dintegrate_dv(),
    Phia_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Phix_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    Phiu_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimu_(robot.dimu()),
    dimi_(0) { 
  if (robot.hasFloatingBase()) {
    dintegrate_dq.resize(robot.dimv(), robot.dimv());
    dintegrate_dv.resize(robot.dimv(), robot.dimv());
    dintegrate_dq.setZero();
    dintegrate_dv.setZero();
  }
}


inline SplitStateConstraintJacobian::SplitStateConstraintJacobian() 
  : dintegrate_dq(),
    dintegrate_dv(),
    Phia_full_(),
    Phix_full_(),
    Phiu_full_(),
    dimv_(0),
    dimx_(0),
    dimu_(0),
    dimi_(0) { 
}


inline SplitStateConstraintJacobian::~SplitStateConstraintJacobian() { 
}


inline void SplitStateConstraintJacobian::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimf();
}


inline void SplitStateConstraintJacobian::setImpulseStatus() {
  dimi_ = 0;
}


inline int SplitStateConstraintJacobian::dimi() const {
  return dimi_;
}


inline Eigen::Block<Eigen::MatrixXd> SplitStateConstraintJacobian::Phia() {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintJacobian::Phia() const {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitStateConstraintJacobian::Phix() {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintJacobian::Phix() const {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitStateConstraintJacobian::Phiq() {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintJacobian::Phiq() const {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitStateConstraintJacobian::Phiv() {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintJacobian::Phiv() const {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SplitStateConstraintJacobian::Phiu() {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintJacobian::Phiu() const {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline bool SplitStateConstraintJacobian::isApprox(
    const SplitStateConstraintJacobian& other) const {
  if (dimi() != other.dimi()) return false;
  if (!Phia().isApprox(other.Phia())) return false;
  if (!Phix().isApprox(other.Phix())) return false;
  if (!Phiu().isApprox(other.Phiu())) return false;
  if (!dintegrate_dq.isApprox(other.dintegrate_dq)) return false;
  if (!dintegrate_dv.isApprox(other.dintegrate_dv)) return false;
  return true;
}


inline bool SplitStateConstraintJacobian::hasNaN() const {
  if (Phia().hasNaN()) return true;
  if (Phix().hasNaN()) return true;
  if (Phiu().hasNaN()) return true;
  if (dintegrate_dq.hasNaN()) return true;
  if (dintegrate_dv.hasNaN()) return true;
  return false;
}

} // namespace idocp

#endif // IDOCP_SPLIT_STATE_CONSTRAINT_JACOBIAN_HXX_