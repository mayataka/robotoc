#ifndef IDOCP_RICCATI_RECURSION_HXX_ 
#define IDOCP_RICCATI_RECURSION_HXX_

#include "idocp/ocp/riccati_recursion.hpp"

#include <cassert>

namespace idocp {

inline RiccatiRecursion::RiccatiRecursion(const Robot& robot, const double T, 
                                          const int N, 
                                          const int max_num_impulse,
                                          const int nproc)
  : N_(N),
    max_num_impulse_(max_num_impulse),
    nproc_(nproc),
    dtau_(T/N),
    riccati_factorizer_(N+1, RiccatiFactorizer(robot),
                        max_num_impulse, ImpulseRiccatiFactorizer(robot)) {
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
    dtau_(0),
    riccati_factorizer_() {
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
    const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual, 
    HybridRiccatiFactorization& riccati_factorization) {
  for (int i=N_-1; i>=0; --i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_impulse;
      assert(dtau_impulse > 0);
      assert(dtau_impulse < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer_.aux[impulse_index].backwardRiccatiRecursion(
          riccati_factorization[i+1], dtau_aux, 
          kkt_matrix.aux[impulse_index],
          kkt_residual.aux[impulse_index],
          riccati_factorization.aux[impulse_index]);
      riccati_factorizer_.impulse[impulse_index].backwardRiccatiRecursion(
          riccati_factorization.aux[impulse_index],  
          kkt_matrix.impulse[impulse_index],
          kkt_residual.impulse[impulse_index],
          riccati_factorization.impulse[impulse_index]);
      riccati_factorizer_[i].backwardRiccatiRecursion(
          riccati_factorization.impulse[impulse_index], dtau_impulse, 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
    else if (contact_sequence.existLiftStage(i)) {
      const int lift_index = contact_sequence.liftIndex(i);
      const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_lift;
      assert(dtau_lift > 0);
      assert(dtau_lift < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer_.lift[lift_index].backwardRiccatiRecursion(
          riccati_factorization[i+1], dtau_aux, 
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
          riccati_factorization.lift[lift_index]);
      riccati_factorizer_[i].backwardRiccatiRecursion(
          riccati_factorization.lift[lift_index], dtau_lift, 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
    else {
      riccati_factorizer_[i].backwardRiccatiRecursion(
          riccati_factorization[i+1], dtau_, 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
  }
}


inline void RiccatiRecursion::forwardRiccatiRecursionParallel(
    const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
    HybridKKTResidual& kkt_residual) {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + N_impulse + N_lift;
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      riccati_factorizer_[i].forwardRiccatiRecursionParallel(
          kkt_matrix[i], kkt_residual[i], exist_state_constraint);
    } 
    else if (i < N_+N_impulse) {
      const int impulse_index = contact_sequence.impulseIndex(i-N_);
      riccati_factorizer_.aux[impulse_index].forwardRiccatiRecursionParallel(
          kkt_matrix.aux[impulse_index], kkt_residual.aux[impulse_index], 
          exist_state_constraint);
    }
    else {
      const int lift_index = contact_sequence.impulseIndex(i-N_-N_impulse);
      riccati_factorizer_.lift[lift_index].forwardRiccatiRecursionParallel(
          kkt_matrix.lift[lift_index], kkt_residual.lift[lift_index], 
          exist_state_constraint);
    }
  }
}


inline void RiccatiRecursion::forwardRiccatiRecursionSerial(
    const ContactSequence& contact_sequence, const HybridKKTMatrix& kkt_matrix, 
    const HybridKKTResidual& kkt_residual, 
    HybridRiccatiFactorization& riccati_factorization) {
  const bool exist_state_constraint = contact_sequence.existImpulseStage();
  for (int i=0; i<N_; ++i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_impulse;
      assert(dtau_impulse > 0);
      assert(dtau_impulse < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      riccati_factorizer_[i].forwardRiccatiRecursionSerial(
          riccati_factorization[i], kkt_matrix[i], kkt_residual[i],
          dtau_impulse, riccati_factorization.impulse[impulse_index],
          exist_state_constraint);
      riccati_factorizer_.impulse[impulse_index].forwardRiccatiRecursionSerial(
          riccati_factorization.impulse[impulse_index],
          kkt_matrix.impulse[impulse_index],
          kkt_residual.impulse[impulse_index],
          riccati_factorization.aux[impulse_index]);
      riccati_factorizer_.aux[impulse_index].forwardRiccatiRecursionSerial(
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
      riccati_factorizer_[i].forwardRiccatiRecursionSerial(
          riccati_factorization[i], kkt_matrix[i], kkt_residual[i], dtau_lift, 
          riccati_factorization.lift[lift_index]);
      riccati_factorizer_.lift[lift_index].forwardRiccatiRecursionSerial(
          riccati_factorization.lift[lift_index], kkt_matrix.lift[lift_index], 
          kkt_residual.lift[lift_index], dtau_aux, riccati_factorization[i+1]);
    }
    else {
      riccati_factorizer_[i].forwardRiccatiRecursionSerial(
        riccati_factorization[i], kkt_matrix[i], kkt_residual[i], dtau_, 
        riccati_factorization[i+1]);
    }
  }
}


inline void RiccatiRecursion::backwardStateConstraintFactorization(
    const ContactSequence& contact_sequence, 
    const HybridKKTMatrix& kkt_matrix, 
    std::vector<StateConstraintRiccatiFactorization>& constraint_factorization) const {
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<num_impulse; ++i) {
    const int impulse_stage = contact_sequence.timeStageBeforeImpulse(i);
    for (int j=impulse_stage; j>=0; --j) {
      if (contact_sequence.existImpulseStage(i)) {
        const int impulse_index = contact_sequence.impulseIndex(i);
        const double dtau_impulse = contact_sequence.impulseTime(impulse_index) - i * dtau_;
        const double dtau_aux = dtau_ - dtau_impulse;
        assert(dtau_impulse > 0);
        assert(dtau_impulse < dtau_);
        assert(dtau_aux > 0);
        assert(dtau_aux < dtau_);
        riccati_factorizer_.aux[impulse_index].backwardStateConstraintFactorization(
            constraint_factorization[i].T(j+1), 
            kkt_matrix.aux[impulse_index], dtau_aux, 
            constraint_factorization[i].T_aux(impulse_index));
        riccati_factorizer_.impulse[impulse_index].backwardStateConstraintFactorization(
            constraint_factorization[i].T_aux(impulse_index), 
            kkt_matrix.impulse[impulse_index], 
            constraint_factorization[i].T_impulse(impulse_index));
        riccati_factorizer_[i].backwardStateConstraintFactorization(
            constraint_factorization[i].T_impulse(impulse_index), 
            kkt_matrix[i], dtau_impulse, 
            constraint_factorization[i].T(j));
      }
      else if (contact_sequence.existLiftStage(i)) {
        const int lift_index = contact_sequence.liftIndex(i);
        const double dtau_lift = contact_sequence.liftTime(lift_index) - i * dtau_;
        const double dtau_aux = dtau_ - dtau_lift;
        assert(dtau_lift > 0);
        assert(dtau_lift < dtau_);
        assert(dtau_aux > 0);
        assert(dtau_aux < dtau_);
        riccati_factorizer_.lift[lift_index].backwardStateConstraintFactorization(
            constraint_factorization[i].T(j+1), 
            kkt_matrix.lift[lift_index], dtau_aux, 
            constraint_factorization[i].T_lift(lift_index));
        riccati_factorizer_[i].backwardStateConstraintFactorization(
            constraint_factorization[i].T_lift(lift_index), 
            kkt_matrix[i], dtau_lift, 
            constraint_factorization[i].T(j));
      }
      else {
        riccati_factorizer_[i].backwardStateConstraintFactorization(
            constraint_factorization[i].T(j+1), 
            kkt_matrix[j], dtau_, 
            constraint_factorization[i].T(j));
      }
    }
  }
}

} // namespace idocp

#endif // IDOCP_RICCATI_RECURSION_HXX_ 