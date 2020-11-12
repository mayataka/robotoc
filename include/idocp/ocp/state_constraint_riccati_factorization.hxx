#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_

#include "idocp/ocp/state_constraint_riccati_factorization.hpp"

#include <cassert>

namespace idocp {

inline StateConstraintRiccatiFactorization::StateConstraintRiccatiFactorization(
    const Robot& robot, const int N, const int max_num_impulse) 
  : T_full_(N, Eigen::MatrixXd(2*robot.dimv(), robot.max_dimf())),
    T_impulse_full_(max_num_impulse, Eigen::MatrixXd(2*robot.dimv(), 
                                                     robot.max_dimf())),
    T_lift_full_(max_num_impulse, Eigen::MatrixXd(2*robot.dimv(), 
                                                  robot.max_dimf())),
    E_full_(Eigen::MatrixXd(robot.max_dimf(), 2*robot.dimv())),
    EN_full_(Eigen::MatrixXd(robot.max_dimf(), 2*robot.dimv())),
    ENEt_full_(Eigen::MatrixXd(robot.max_dimf(), robot.max_dimf())),
    e_full_(Eigen::VectorXd(robot.max_dimf())),
    N_(N),
    max_num_impulse_(max_num_impulse),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimf_(0) { 
}


inline StateConstraintRiccatiFactorization::
StateConstraintRiccatiFactorization() 
  : T_full_(),
    T_impulse_full_(),
    T_lift_full_(),
    E_full_(),
    EN_full_(),
    ENEt_full_(),
    e_full_(),
    N_(0),
    max_num_impulse_(0),
    dimv_(0),
    dimx_(0),
    dimf_(0) { 
}


inline StateConstraintRiccatiFactorization::
~StateConstraintRiccatiFactorization() { 
}


inline void StateConstraintRiccatiFactorization::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimp();
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T(const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return T_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return T_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_impulse(const int impulse_index) {
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return T_impulse_full_[impulse_index].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_impulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return T_impulse_full_[impulse_index].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_lift(const int lift_index) {
  assert(lift_index >= 0);
  assert(lift_index < max_num_impulse_);
  return T_lift_full_[lift_index].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < max_num_impulse_);
  return T_lift_full_[lift_index].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::Eq() {
  return E_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::Eq() const {
  return E_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::EN() {
  return EN_full_.topLeftCorner(dimf_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::EN() const {
  return EN_full_.topLeftCorner(dimf_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENq() {
  return EN_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENq() const {
  return EN_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENEt() {
  return ENEt_full_.topLeftCorner(dimf_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENEt() const {
  return ENEt_full_.topLeftCorner(dimf_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
StateConstraintRiccatiFactorization::e() {
  return e_full_.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
StateConstraintRiccatiFactorization::e() const {
  return e_full_.head(dimf_);
}
} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_ 