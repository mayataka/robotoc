#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_

#include "idocp/ocp/state_constraint_riccati_factorization.hpp"

#include <cassert>

namespace idocp {

inline StateConstraintRiccatiFactorization::StateConstraintRiccatiFactorization(
    const Robot& robot, const int N) 
  : T_full_(N, Eigen::MatrixXd(2*robot.dimv(), robot.max_dimf())),
    T_aux_full_(N, Eigen::MatrixXd(2*robot.dimv(), robot.max_dimf())),
    EqNqq_full_(Eigen::MatrixXd(robot.max_dimf(), robot.dimv())),
    ENEt_full_(Eigen::MatrixXd(robot.max_dimf(), robot.max_dimf())),
    N_(N),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimf_(0),
    is_active_(false) { 
}


inline StateConstraintRiccatiFactorization::
StateConstraintRiccatiFactorization() 
  : T_full_(),
    T_aux_full_(),
    EqNqq_full_(),
    ENEt_full_(),
    N_(0),
    dimv_(0),
    dimx_(0),
    dimf_(0), 
    is_active_(false) { 
}


inline StateConstraintRiccatiFactorization::
~StateConstraintRiccatiFactorization() { 
}


inline void StateConstraintRiccatiFactorization::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimp();
  is_active_ = impulse_status.hasActiveImpulse();
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T(const int time_stage) {
  return T_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T(const int time_stage) const {
  return T_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_aux(const int time_stage) {
  return T_aux_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_aux(const int time_stage) const {
  return T_aux_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENEt() {
  return ENEt_full_.topLeftCorner(dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENEt() const {
  return ENEt_full_.topLeftCorner(dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::EqNqq() {
  return EqNqq_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::EqNqq() const {
  return EqNqq_full_.topLeftCorner(dimf_, dimv_);
}

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_ 