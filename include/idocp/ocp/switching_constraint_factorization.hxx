#ifndef IDOCP_SWITCHING_CONSTRAINT_FACTORIZATION_HXX_ 
#define IDOCP_SWITCHING_CONSTRAINT_FACTORIZATION_HXX_

#include "idocp/ocp/switching_constraint_factorization.hpp"

#include <cassert>

namespace idocp {

inline SwitchingConstraintFactorization::
SwitchingConstraintFactorization(const Robot& robot) 
  : q(Eigen::VectorXd::Zero(robot.dimq())),
    dq(Eigen::VectorXd::Zero(robot.dimv())),
    dintegrate_dq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    dintegrate_dv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Phia_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    Phix_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    Phiu_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    DGinv_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    S_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Sinv_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    SinvDGinv_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimu())),
    M_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    m_full_(Eigen::VectorXd::Zero(robot.max_dimf())),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimu_(robot.dimu()),
    dimi_(0) { 
}


inline SwitchingConstraintFactorization::
SwitchingConstraintFactorization() 
  : q(),
    dq(),
    dintegrate_dq(),
    dintegrate_dv(),
    Phia_full_(),
    Phix_full_(),
    Phiu_full_(),
    DGinv_full_(),
    S_full_(),
    Sinv_full_(),
    SinvDGinv_full_(),
    M_full_(),
    m_full_(),
    dimv_(0),
    dimx_(0),
    dimu_(0),
    dimi_(0) { 
}


inline SwitchingConstraintFactorization::~SwitchingConstraintFactorization() { 
}


inline void SwitchingConstraintFactorization::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimi_ = impulse_status.dimf();
}


inline void SwitchingConstraintFactorization::setImpulseStatus() {
  dimi_ = 0;
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::Phia() {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::Phia() const {
  return Phia_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::Phix() {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::Phix() const {
  return Phix_full_.topLeftCorner(dimi_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::Phiq() {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::Phiq() const {
  return Phix_full_.topLeftCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::Phiv() {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::Phiv() const {
  return Phix_full_.topRightCorner(dimi_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::Phiu() {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::Phiu() const {
  return Phiu_full_.topLeftCorner(dimi_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::DGinv() {
  return DGinv_full_.topLeftCorner(dimi_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::DGinv() const {
  return DGinv_full_.topLeftCorner(dimi_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::S() {
  return S_full_.topLeftCorner(dimi_, dimi_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::S() const {
  return S_full_.topLeftCorner(dimi_, dimi_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::Sinv() {
  return Sinv_full_.topLeftCorner(dimi_, dimi_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::Sinv() const {
  return Sinv_full_.topLeftCorner(dimi_, dimi_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SwitchingConstraintFactorization::SinvDGinv() {
  return SinvDGinv_full_.topLeftCorner(dimi_, dimu_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::SinvDGinv() const {
  return SinvDGinv_full_.topLeftCorner(dimi_, dimu_);
}


inline Eigen::Block<Eigen::MatrixXd> SwitchingConstraintFactorization::M() {
  return M_full_.topLeftCorner(dimi_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SwitchingConstraintFactorization::M() const {
  return M_full_.topLeftCorner(dimi_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SwitchingConstraintFactorization::m() {
  return m_full_.head(dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SwitchingConstraintFactorization::m() const {
  return m_full_.head(dimi_);
}

} // namespace idocp

#endif // IDOCP_SWITCHING_CONSTRAINT_FACTORIZATION_HXX_