#ifndef IDOCP_RICCATI_RECURSION_HXX_
#define IDOCP_RICCATI_RECURSION_HXX_

#include "idocp/ocp/riccati_recursion.hpp"

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
    dtau_(T/N),
    riccati_factorization_(N+1, RiccatiFactorization(robot)),
    impulse_riccati_factorization_(max_num_impulse, RiccatiFactorization(robot)),
    aux_riccati_factorization_(max_num_impulse, RiccatiFactorization(robot)),
    lift_riccati_factorization_(max_num_impulse, RiccatiFactorization(robot)),
    constraint_factorization_(
        max_num_impulse, 
        StateConstraintRiccatiFactorization(robot, N, max_num_impulse)),
    constraint_factorizer_(robot, max_num_impulse, nproc),
    primal_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)), 
    dual_step_sizes_(Eigen::VectorXd::Zero(N+1+3*max_num_impulse)) {
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
    riccati_factorization_(),
    impulse_riccati_factorization_(),
    aux_riccati_factorization_(),
    lift_riccati_factorization_(),
    constraint_factorization_(),
    constraint_factorizer_(),
    primal_step_sizes_(), 
    dual_step_sizes_() {
}


inline RiccatiRecursion::~RiccatiRecursion() {
}


inline void RiccatiRecursion::backwardRiccatiRecursionTerminal(
    const TerminalOCP& terminal_ocp) {
  terminal_ocp.backwardRiccatiRecursion(riccati_factorization_[N_]);
}


inline void RiccatiRecursion::backwardRiccatiRecursion(
    const ContactSequence& contact_sequence, HybridSplitOCPs& split_ocps) {
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
      split_ocps.aux[impulse_index].backwardRiccatiRecursion(
          riccati_factorization_[i+1], dtau_aux, 
          aux_riccati_factorization_[impulse_index]);
      split_ocps.impulse[impulse_index].backwardRiccatiRecursion(
          aux_riccati_factorization_[impulse_index],  
          impulse_riccati_factorization_[impulse_index]);
      split_ocps[i].backwardRiccatiRecursion(
          impulse_riccati_factorization_[impulse_index], dtau_impulse, 
          riccati_factorization_[i]);
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
      split_ocps.lift[lift_index].backwardRiccatiRecursion(
          riccati_factorization_[i+1], dtau_aux, 
          lift_riccati_factorization_[lift_index]);
      split_ocps[i].backwardRiccatiRecursion(
          impulse_riccati_factorization_[lift_index], dtau_lift, 
          riccati_factorization_[i]);
    }
    else {
      split_ocps[i].backwardRiccatiRecursion(riccati_factorization_[i+1], 
                                              dtau_, riccati_factorization_[i]);
    }
  }
}


inline void RiccatiRecursion::forwardRiccatiRecursionParallel(
    const ContactSequence& contact_sequence, 
    HybridSplitOCPs& split_ocps) const {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 2 * N_impulse + N_lift;
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      split_ocps[i].forwardRiccatiRecursionParallel();
    } 
    else if (i < N_+N_impulse) {
      // split_ocps.impulse[i-N_].forwardRiccatiRecursionParallel();
    }
    else if (i < N_+2*N_impulse) {
      split_ocps.aux[i-N_].forwardRiccatiRecursionParallel();
    }
    else {
      split_ocps.lift[i-N_].forwardRiccatiRecursionParallel();
    }
  }
}


inline void RiccatiRecursion::forwardRiccatiRecursionSerial(
    const ContactSequence& contact_sequence, 
    HybridSplitOCPs& split_ocps) {
  for (int i=0; i<N_; --i) {
    if (contact_sequence.existImpulseStage(i)) {
      const int impulse_index = contact_sequence.impulseIndex(i);
      const double dtau_impulse
          = contact_sequence.impulseTime(impulse_index) - i * dtau_;
      const double dtau_aux = dtau_ - dtau_impulse;
      assert(dtau_impulse > 0);
      assert(dtau_impulse < dtau_);
      assert(dtau_aux > 0);
      assert(dtau_aux < dtau_);
      split_ocps[i].forwardRiccatiRecursionSerial(
          riccati_factorization_[i], dtau_impulse, 
          impulse_riccati_factorization_[impulse_index]);
      split_ocps.impulse[impulse_index].forwardRiccatiRecursionSerial(
          impulse_riccati_factorization_[impulse_index],
          aux_riccati_factorization_[impulse_index]);
      split_ocps.aux[impulse_index].forwardRiccatiRecursionSerial(
          aux_riccati_factorization_[impulse_index], dtau_aux, 
          riccati_factorization_[i+1]);
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
      split_ocps[i].forwardRiccatiRecursionSerial(
          riccati_factorization_[i], dtau_lift, 
          lift_riccati_factorization_[lift_index]);
      split_ocps.lift[i].forwardRiccatiRecursionSerial(
          lift_riccati_factorization_[lift_index], dtau_aux, 
          riccati_factorization_[i+1]);
    }
    else {
      split_ocps[i].forwardRiccatiRecursionSerial(riccati_factorization_[i], 
                                                  dtau_, 
                                                  riccati_factorization_[i+1]);
    }
  }
}


template <typename VectorType>
inline void RiccatiRecursion::computeLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence, 
    const HybridSplitOCPs& split_ocps,
    const Eigen::MatrixBase<VectorType>& dx0,
    HybridSplitDirections& d) {
  constraint_factorizer_.computeLagrangeMultiplierDirection(
      contact_sequence, impulse_riccati_factorization_, 
      constraint_factorization_, dx0, d.impulse);
}


inline void RiccatiRecursion::computeDirectionAndStepSize(
    std::vector<Robot>& robots, const ContactSequence& contact_sequence, 
    const HybridSplitOCPs& split_ocps, const HybridSplitSolutions& s, 
    HybridSplitDirections& d) {
  
}


inline double RiccatiRecursion::primalStepSize() const {
  return primal_step_sizes_.minCoeff();
}


inline double RiccatiRecursion::dualStepSize() const {
  return dual_step_sizes_.minCoeff();
}

} // namespace idocp 

#endif // IDOCP_RICCATI_RECURSION_HXX_ 