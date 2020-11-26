#ifndef IDOCP_RICCATI_RECURSION_HXX_ 
#define IDOCP_RICCATI_RECURSION_HXX_

#include "idocp/ocp/riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp {

inline RiccatiRecursion::RiccatiRecursion(const Robot& robot, const double T, 
                                          const int N, 
                                          const int max_num_impulse,
                                          const int nproc)
  : N_(N),
    max_num_impulse_(max_num_impulse),
    nproc_(nproc),
    dimv_(robot.dimv()),
    dtau_(T/N) {
  try {
    if (T <= 0) {
      throw std::out_of_range("invalid value: T must be positive!");
    }
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
    }
    if (max_num_impulse < 0) {
      throw std::out_of_range("invalid value: max_num_impulse must be non-negative!");
    }
    if (nproc <= 0) {
      throw std::out_of_range("invalid value: nproc must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline RiccatiRecursion::RiccatiRecursion()
  : N_(0),
    max_num_impulse_(0),
    nproc_(0),
    dimv_(0),
    dtau_(0) {
}


inline RiccatiRecursion::~RiccatiRecursion() {
}


inline void RiccatiRecursion::backwardRiccatiRecursionTerminal(
    const HybridKKTMatrix& kkt_matrix, const HybridKKTResidual& kkt_residual,
    HybridRiccatiFactorization& riccati_factorization) const {
  riccati_factorization[N_].Pqq = kkt_matrix[N_].Qqq();
  riccati_factorization[N_].Pvv = kkt_matrix[N_].Qvv();
  riccati_factorization[N_].sq = - kkt_residual[N_].lq();
  riccati_factorization[N_].sv = - kkt_residual[N_].lv();
}


inline void RiccatiRecursion::backwardRiccatiRecursion(
    HybridRiccatiFactorizer& riccati_factorizer,
    const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual, 
    HybridRiccatiFactorization& riccati_factorization) {
  for (int i=N_-1; i>=0; --i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse 
          = contact_sequence.impulseTime(impulse_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_impulse;
      assert(dtau_impulse > 0);
      assert(dtau_impulse < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer.aux[impulse_index].backwardRiccatiRecursion(
          riccati_factorization[i+1], dtau_aux, 
          kkt_matrix.aux[impulse_index],
          kkt_residual.aux[impulse_index],
          riccati_factorization.aux[impulse_index]);
      riccati_factorizer.impulse[impulse_index].backwardRiccatiRecursion(
          riccati_factorization.aux[impulse_index],  
          kkt_matrix.impulse[impulse_index],
          kkt_residual.impulse[impulse_index],
          riccati_factorization.impulse[impulse_index]);
      riccati_factorizer[i].backwardRiccatiRecursion(
          riccati_factorization.impulse[impulse_index], dtau_impulse, 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift 
          = contact_sequence.liftTime(lift_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_lift;
      assert(dtau_lift > 0);
      assert(dtau_lift < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer.lift[lift_index].backwardRiccatiRecursion(
          riccati_factorization[i+1], dtau_aux, 
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
          riccati_factorization.lift[lift_index]);
      riccati_factorizer[i].backwardRiccatiRecursion(
          riccati_factorization.lift[lift_index], dtau_lift, 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
    else {
      riccati_factorizer[i].backwardRiccatiRecursion(
          riccati_factorization[i+1], dtau_, 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
  }
}


inline void RiccatiRecursion::forwardRiccatiRecursionParallel(
    HybridRiccatiFactorizer& riccati_factorizer,
    const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual,
      std::vector<StateConstraintRiccatiFactorization>& constraint_factorization) {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 2*N_impulse + N_lift;
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      riccati_factorizer[i].forwardRiccatiRecursionParallel(
          kkt_matrix[i], kkt_residual[i], exist_state_constraint);
    } 
    else if (i < N_+N_impulse) {
      const int impulse_index = i - N_;
      constraint_factorization[impulse_index].Eq() 
          = kkt_matrix.impulse[impulse_index].Pq();
      constraint_factorization[impulse_index].e() 
          = kkt_residual.impulse[impulse_index].P();
      constraint_factorization[impulse_index].T_impulse(impulse_index).topRows(dimv_)
          = kkt_matrix.impulse[impulse_index].Pq().transpose();
      constraint_factorization[impulse_index].T_impulse(impulse_index).bottomRows(dimv_).setZero();
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_index = i - N_ - N_impulse;
      riccati_factorizer.aux[impulse_index].forwardRiccatiRecursionParallel(
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
          exist_state_constraint);
    }
    else {
      const int lift_index = i - N_ - 2 * N_impulse;
      riccati_factorizer.lift[lift_index].forwardRiccatiRecursionParallel(
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
          exist_state_constraint);
    }
  }
}


inline void RiccatiRecursion::forwardStateConstraintFactorization(
    HybridRiccatiFactorizer& riccati_factorizer,
    const ContactSequence& contact_sequence, const HybridKKTMatrix& kkt_matrix, 
    const HybridKKTResidual& kkt_residual, 
    HybridRiccatiFactorization& riccati_factorization) {
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  riccati_factorizer[0].forwardStateConstraintFactorizationInitial(
      riccati_factorization[0]);
  for (int i=0; i<N_; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse 
          = contact_sequence.impulseTime(impulse_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_impulse;
      assert(dtau_impulse > 0);
      assert(dtau_impulse < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer[i].forwardStateConstraintFactorization(
          riccati_factorization[i], kkt_matrix[i], kkt_residual[i],
          dtau_impulse, riccati_factorization.impulse[impulse_index],
          exist_state_constraint);
      riccati_factorizer.impulse[impulse_index].forwardStateConstraintFactorization(
          riccati_factorization.impulse[impulse_index],
          kkt_matrix.impulse[impulse_index],
          kkt_residual.impulse[impulse_index],
          riccati_factorization.aux[impulse_index], exist_state_constraint);
      riccati_factorizer.aux[impulse_index].forwardStateConstraintFactorization(
          riccati_factorization.aux[impulse_index],
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index],
          dtau_aux, riccati_factorization[i+1], exist_state_constraint);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_lift;
      assert(dtau_lift > 0);
      assert(dtau_lift < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer[i].forwardStateConstraintFactorization(
          riccati_factorization[i], kkt_matrix[i], kkt_residual[i], dtau_lift, 
          riccati_factorization.lift[lift_index], exist_state_constraint);
      riccati_factorizer.lift[lift_index].forwardStateConstraintFactorization(
          riccati_factorization.lift[lift_index], kkt_matrix.lift[lift_index], 
          kkt_residual.lift[lift_index], dtau_aux, riccati_factorization[i+1],
          exist_state_constraint);
    }
    else {
      riccati_factorizer[i].forwardStateConstraintFactorization(
        riccati_factorization[i], kkt_matrix[i], kkt_residual[i], dtau_, 
        riccati_factorization[i+1], exist_state_constraint);
    }
  }
}


inline void RiccatiRecursion::backwardStateConstraintFactorization(
    const HybridRiccatiFactorizer& riccati_factorizer,
    const ContactSequence& contact_sequence, const HybridKKTMatrix& kkt_matrix, 
    std::vector<StateConstraintRiccatiFactorization>& constraint_factorization) const {
  const int num_constraint = contact_sequence.totalNumImpulseStages();
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<num_constraint; ++i) {
    backwardStateConstraintFactorization(riccati_factorizer, contact_sequence, 
                                         kkt_matrix, 
                                         constraint_factorization[i], i);
  }
}


inline void RiccatiRecursion::backwardStateConstraintFactorization(
    const HybridRiccatiFactorizer& riccati_factorizer,
    const ContactSequence& contact_sequence, 
    const HybridKKTMatrix& kkt_matrix, 
    StateConstraintRiccatiFactorization& constraint_factorization,
    const int constraint_index) const {
  assert(constraint_index >= 0);
  assert(constraint_index < contact_sequence.totalNumImpulseStages());
  const int time_stage_before_constraint 
      = contact_sequence.timeStageBeforeImpulse(constraint_index);
  const double dtau_constraint 
      = contact_sequence.impulseTime(constraint_index) 
          - time_stage_before_constraint * dtau_;
  riccati_factorizer[time_stage_before_constraint].backwardStateConstraintFactorization(
      constraint_factorization.T_impulse(constraint_index), 
      kkt_matrix[time_stage_before_constraint], dtau_constraint, 
      constraint_factorization.T(time_stage_before_constraint));
  for (int i=time_stage_before_constraint-1; i>=0; --i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse 
          = contact_sequence.impulseTime(impulse_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_impulse;
      assert(dtau_impulse > 0);
      assert(dtau_impulse < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer.aux[impulse_index].backwardStateConstraintFactorization(
          constraint_factorization.T(i+1), 
          kkt_matrix.aux[impulse_index], dtau_aux, 
          constraint_factorization.T_aux(impulse_index));
      riccati_factorizer.impulse[impulse_index].backwardStateConstraintFactorization(
          constraint_factorization.T_aux(impulse_index), 
          kkt_matrix.impulse[impulse_index], 
          constraint_factorization.T_impulse(impulse_index));
      riccati_factorizer[i].backwardStateConstraintFactorization(
          constraint_factorization.T_impulse(impulse_index), 
          kkt_matrix[i], dtau_impulse, 
          constraint_factorization.T(i));
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift 
          = contact_sequence.liftTime(lift_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_lift;
      assert(dtau_lift > 0);
      assert(dtau_lift < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer.lift[lift_index].backwardStateConstraintFactorization(
          constraint_factorization.T(i+1), 
          kkt_matrix.lift[lift_index], dtau_aux, 
          constraint_factorization.T_lift(lift_index));
      riccati_factorizer[i].backwardStateConstraintFactorization(
          constraint_factorization.T_lift(lift_index), 
          kkt_matrix[i], dtau_lift, 
          constraint_factorization.T(i));
    }
    else {
      riccati_factorizer[i].backwardStateConstraintFactorization(
          constraint_factorization.T(i+1), 
          kkt_matrix[i], dtau_, 
          constraint_factorization.T(i));
    }
  }
}


inline void RiccatiRecursion::forwardRiccatiRecursion(
    const HybridRiccatiFactorizer& riccati_factorizer,
    const ContactSequence& contact_sequence, const HybridKKTMatrix& kkt_matrix, 
    const HybridKKTResidual& kkt_residual, 
    const HybridRiccatiFactorization& riccati_factorization, 
    HybridDirection& d) {
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  for (int i=0; i<N_; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse 
          = contact_sequence.impulseTime(impulse_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_impulse;
      assert(dtau_impulse > 0);
      assert(dtau_impulse < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer[i].forwardRiccatiRecursion(
          kkt_matrix[i], kkt_residual[i],
          riccati_factorization.impulse[impulse_index], 
          d[i], dtau_impulse, d.impulse[impulse_index],
          exist_state_constraint);
      riccati_factorizer.impulse[impulse_index].forwardRiccatiRecursion(
          kkt_matrix.impulse[impulse_index], 
          kkt_residual.impulse[impulse_index],
          d.impulse[impulse_index], d.aux[impulse_index]);
      riccati_factorizer.aux[impulse_index].forwardRiccatiRecursion(
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index],
          riccati_factorization[i+1], d.aux[impulse_index], dtau_aux, d[i+1], 
          exist_state_constraint);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift 
          = contact_sequence.liftTime(lift_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_lift;
      assert(dtau_lift > 0);
      assert(dtau_lift < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer[i].forwardRiccatiRecursion(
          kkt_matrix[i], kkt_residual[i], riccati_factorization.lift[lift_index], 
          d[i], dtau_lift, d.lift[lift_index], exist_state_constraint);
      riccati_factorizer.lift[lift_index].forwardRiccatiRecursion(
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index],
          riccati_factorization[i+1], d.lift[lift_index], dtau_aux, d[i+1], 
          exist_state_constraint);
    }
    else {
      riccati_factorizer[i].forwardRiccatiRecursion(
        kkt_matrix[i], kkt_residual[i], riccati_factorization[i+1], 
        d[i], dtau_, d[i+1], exist_state_constraint);
    }
  }
}

} // namespace idocp

#endif // IDOCP_RICCATI_RECURSION_HXX_ 