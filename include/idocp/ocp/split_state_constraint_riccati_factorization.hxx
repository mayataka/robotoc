#ifndef IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_
#define IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_

#include "idocp/ocp/split_state_constraint_riccati_factorization.hpp"

#include <stdexcept>
#include <cassert>

namespace idocp {

inline SplitStateConstraintRiccatiFactorization::
SplitStateConstraintRiccatiFactorization(const Robot& robot, const int N, 
                                         const int max_num_impulse) 
  : T_full_(N, Eigen::MatrixXd(2*robot.dimv(), robot.max_dimf())),
    T_impulse_full_(max_num_impulse, Eigen::MatrixXd(2*robot.dimv(), 
                                                     robot.max_dimf())),
    T_aux_full_(max_num_impulse, Eigen::MatrixXd(2*robot.dimv(), 
                                                 robot.max_dimf())),
    T_lift_full_(max_num_impulse, Eigen::MatrixXd(2*robot.dimv(), 
                                                  robot.max_dimf())),
    E_full_(Eigen::MatrixXd(robot.max_dimf(), 2*robot.dimv())),
    EN_full_(Eigen::MatrixXd(robot.max_dimf(), 2*robot.dimv())),
    N_(N),
    max_num_impulse_(max_num_impulse),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimf_(0) { 
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline SplitStateConstraintRiccatiFactorization::
SplitStateConstraintRiccatiFactorization() 
  : T_full_(),
    T_impulse_full_(),
    T_aux_full_(),
    T_lift_full_(),
    E_full_(),
    EN_full_(),
    N_(0),
    max_num_impulse_(0),
    dimv_(0),
    dimx_(0),
    dimf_(0) { 
}


inline SplitStateConstraintRiccatiFactorization::
~SplitStateConstraintRiccatiFactorization() { 
}


inline void SplitStateConstraintRiccatiFactorization::setImpulseStatus(
    const ImpulseStatus& impulse_status) {
  dimf_ = impulse_status.dimp();
}


inline void SplitStateConstraintRiccatiFactorization::resetImpulseStatus() {
  dimf_ = 0;
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::T(const int time_stage) {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return T_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::T(const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return T_full_[time_stage].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::T_impulse(const int impulse_index) {
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return T_impulse_full_[impulse_index].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::
T_impulse(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return T_impulse_full_[impulse_index].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::T_aux(const int impulse_index) {
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return T_aux_full_[impulse_index].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::T_aux(const int impulse_index) const {
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return T_aux_full_[impulse_index].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::T_lift(const int lift_index) {
  assert(lift_index >= 0);
  assert(lift_index < max_num_impulse_);
  return T_lift_full_[lift_index].topLeftCorner(dimx_, dimf_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::T_lift(const int lift_index) const {
  assert(lift_index >= 0);
  assert(lift_index < max_num_impulse_);
  return T_lift_full_[lift_index].topLeftCorner(dimx_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::Eq() {
  return E_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::Eq() const {
  return E_full_.topLeftCorner(dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::EN() {
  return EN_full_.topLeftCorner(dimf_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::EN() const {
  return EN_full_.topLeftCorner(dimf_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::ENq() {
  return EN_full_.topLeftCorner(dimf_, dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitStateConstraintRiccatiFactorization::ENq() const {
  return EN_full_.topLeftCorner(dimf_, dimv_);
}

} // namespace idocp

#endif // IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_ 