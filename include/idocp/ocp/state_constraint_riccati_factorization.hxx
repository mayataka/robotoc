#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_

#include "idocp/ocp/state_constraint_riccati_factorization.hpp"

#include <cassert>

namespace idocp {

inline StateConstraintRiccatiFactorization::StateConstraintRiccatiFactorization(
    const Robot& robot, const int N, const int max_num_impulse) 
  : ENT_full_(Eigen::MatrixXd::Zero(max_num_impulse*robot.max_dimf(), 
                                    max_num_impulse*robot.max_dimf())),
    e_full_(Eigen::VectorXd::Zero(max_num_impulse*robot.max_dimf())),
    dxi_full_(Eigen::VectorXd::Zero(max_num_impulse*robot.max_dimf())),
    factorizations_(
        max_num_impulse, 
        SplitStateConstraintRiccatiFactorization(robot, N, max_num_impulse)),
    dimf_(max_num_impulse, 0),
    f_begin_(max_num_impulse, 0),
    N_(N),
    max_num_impulse_(max_num_impulse),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimf_total_(0) { 
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


inline StateConstraintRiccatiFactorization::
StateConstraintRiccatiFactorization() 
  : ENT_full_(),
    e_full_(),
    dxi_full_(),
    factorizations_(),
    dimf_(),
    f_begin_(),
    N_(0),
    max_num_impulse_(0),
    dimv_(0),
    dimx_(0),
    dimf_total_(0) { 
}


inline StateConstraintRiccatiFactorization::
~StateConstraintRiccatiFactorization() { 
}


inline void StateConstraintRiccatiFactorization::setConstraintStatus(
    const ContactSequence& contact_sequence) {
  const int total_num_impulse = contact_sequence.totalNumImpulseStages();
  for (int i=0; i<total_num_impulse; ++i) {
    dimf_[i] = contact_sequence.impulseStatus(i).dimp();
    factorizations_[i].setImpulseStatus(contact_sequence.impulseStatus(i));
  }
  for (int i=total_num_impulse; i<max_num_impulse_; ++i) {
    dimf_[i] = 0;
    factorizations_[i].resetImpulseStatus();
  }
  dimf_total_ = 0;
  for (int i=0; i<total_num_impulse; ++i) {
    f_begin_[i] = dimf_total_;
    dimf_total_ += dimf_[i];
  }
  for (int i=total_num_impulse; i<max_num_impulse_; ++i) {
    f_begin_[i] = 0;
  }
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T(const int constraint_index, 
                                       const int time_stage) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return factorizations_[constraint_index].T(time_stage);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T(const int constraint_index,
                                       const int time_stage) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return factorizations_[constraint_index].T(time_stage);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_impulse(const int constraint_index,
                                               const int impulse_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return factorizations_[constraint_index].T_impulse(impulse_index);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_impulse(const int constraint_index,
                                               const int impulse_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return factorizations_[constraint_index].T_impulse(impulse_index);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_aux(const int constraint_index,
                                           const int impulse_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return factorizations_[constraint_index].T_aux(impulse_index);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_aux(const int constraint_index, 
                                           const int impulse_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(impulse_index >= 0);
  assert(impulse_index < max_num_impulse_);
  return factorizations_[constraint_index].T_aux(impulse_index);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_lift(const int constraint_index,
                                            const int lift_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(lift_index >= 0);
  assert(lift_index < max_num_impulse_);
  return factorizations_[constraint_index].T_lift(lift_index);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::T_lift(const int constraint_index, 
                                            const int lift_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  assert(lift_index >= 0);
  assert(lift_index < max_num_impulse_);
  return factorizations_[constraint_index].T_lift(lift_index);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::Eq(const int constraint_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return factorizations_[constraint_index].Eq();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::Eq(const int constraint_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return factorizations_[constraint_index].Eq();
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::EN(const int constraint_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return factorizations_[constraint_index].EN();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::EN(const int constraint_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return factorizations_[constraint_index].EN();
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENq(const int constraint_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return factorizations_[constraint_index].ENq();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENq(const int constraint_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return factorizations_[constraint_index].ENq();
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENEt(const int constraint_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return ENT(constraint_index, constraint_index);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENEt(const int constraint_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return ENT(constraint_index, constraint_index);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENT(const int constraint_index, 
                                         const int impulse_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return ENT_full_.block(f_begin_[constraint_index], f_begin_[impulse_index], 
                         dimf_[constraint_index], dimf_[impulse_index]);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENT(const int constraint_index, 
                                         const int impulse_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return ENT_full_.block(f_begin_[constraint_index], f_begin_[impulse_index], 
                         dimf_[constraint_index], dimf_[impulse_index]);
}


inline Eigen::Block<Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENT() {
  return ENT_full_.topLeftCorner(dimf_total_, dimf_total_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
StateConstraintRiccatiFactorization::ENT() const {
  return ENT_full_.topLeftCorner(dimf_total_, dimf_total_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
StateConstraintRiccatiFactorization::e(const int constraint_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return e_full_.segment(f_begin_[constraint_index], dimf_[constraint_index]);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
StateConstraintRiccatiFactorization::e(const int constraint_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return e_full_.segment(f_begin_[constraint_index], dimf_[constraint_index]);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
StateConstraintRiccatiFactorization::e() {
  return e_full_.head(dimf_total_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
StateConstraintRiccatiFactorization::e() const {
  return e_full_.head(dimf_total_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
StateConstraintRiccatiFactorization::dxi() {
  return dxi_full_.head(dimf_total_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
StateConstraintRiccatiFactorization::dxi() const {
  return dxi_full_.head(dimf_total_);
}

inline Eigen::VectorBlock<Eigen::VectorXd> 
StateConstraintRiccatiFactorization::dxi(const int constraint_index) {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return dxi_full_.segment(f_begin_[constraint_index], dimf_[constraint_index]);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
StateConstraintRiccatiFactorization::dxi(const int constraint_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < max_num_impulse_);
  return dxi_full_.segment(f_begin_[constraint_index], dimf_[constraint_index]);
}

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HXX_ 