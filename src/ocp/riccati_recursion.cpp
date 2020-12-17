#include "idocp/ocp/riccati_recursion.hpp"

#include <omp.h>
#include <stdexcept>
#include <cassert>

namespace idocp {

RiccatiRecursion::RiccatiRecursion(const Robot& robot, const int N, 
                                   const int nproc)
  : N_(N),
    nproc_(nproc),
    dimv_(robot.dimv()) {
  try {
    if (N <= 0) {
      throw std::out_of_range("invalid value: N must be positive!");
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


RiccatiRecursion::RiccatiRecursion()
  : N_(0),
    nproc_(0),
    dimv_(0) {
}


RiccatiRecursion::~RiccatiRecursion() {
}


void RiccatiRecursion::backwardRiccatiRecursionTerminal(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual,
    RiccatiFactorization& riccati_factorization) const {
  riccati_factorization[N_].Pqq = kkt_matrix[N_].Qqq();
  riccati_factorization[N_].Pvv = kkt_matrix[N_].Qvv();
  riccati_factorization[N_].sq = - kkt_residual[N_].lq();
  riccati_factorization[N_].sv = - kkt_residual[N_].lv();
}


void RiccatiRecursion::backwardRiccatiRecursion(
    RiccatiFactorizer& riccati_factorizer, 
    const OCPDiscretizer& ocp_discretizer, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual, RiccatiFactorization& riccati_factorization) {
  for (int i=N_-1; i>=0; --i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_idx = ocp_discretizer.impulseIndex(i);
      riccati_factorizer.aux[impulse_idx].backwardRiccatiRecursion(
          riccati_factorization[i+1], ocp_discretizer.dtau_aux(impulse_idx), 
          kkt_matrix.aux[impulse_idx], kkt_residual.aux[impulse_idx], 
          riccati_factorization.aux[impulse_idx]);
      riccati_factorizer.impulse[impulse_idx].backwardRiccatiRecursion(
          riccati_factorization.aux[impulse_idx],  
          kkt_matrix.impulse[impulse_idx], kkt_residual.impulse[impulse_idx],
          riccati_factorization.impulse[impulse_idx]);
      riccati_factorizer[i].backwardRiccatiRecursion(
          riccati_factorization.impulse[impulse_idx], ocp_discretizer.dtau(i), 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_idx = ocp_discretizer.liftIndex(i);
      riccati_factorizer.lift[lift_idx].backwardRiccatiRecursion(
          riccati_factorization[i+1], ocp_discretizer.dtau_lift(lift_idx), 
          kkt_matrix.lift[lift_idx], kkt_residual.lift[lift_idx], 
          riccati_factorization.lift[lift_idx]);
      riccati_factorizer[i].backwardRiccatiRecursion(
          riccati_factorization.lift[lift_idx], ocp_discretizer.dtau(i), 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
    else {
      riccati_factorizer[i].backwardRiccatiRecursion(
          riccati_factorization[i+1], ocp_discretizer.dtau(i), 
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i]);
    }
  }
}


void RiccatiRecursion::forwardRiccatiRecursionParallel(
    RiccatiFactorizer& riccati_factorizer, 
    const OCPDiscretizer& ocp_discretizer, KKTMatrix& kkt_matrix, 
    KKTResidual& kkt_residual, 
    StateConstraintRiccatiFactorization& constraint_factorization) {
  const int N_impulse = ocp_discretizer.numImpulseStages();
  const int N_lift = ocp_discretizer.numLiftStages();
  const int N_all = N_ + 2*N_impulse + N_lift;
  const bool exist_state_constraint = ocp_discretizer.existImpulse();
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      riccati_factorizer[i].forwardRiccatiRecursionParallel(
          kkt_matrix[i], kkt_residual[i], exist_state_constraint);
    } 
    else if (i < N_+N_impulse) {
      const int impulse_idx = i - N_;
      constraint_factorization.Eq(impulse_idx) 
          = kkt_matrix.impulse[impulse_idx].Pq();
      constraint_factorization.e(impulse_idx) 
          = kkt_residual.impulse[impulse_idx].P();
      constraint_factorization.T_impulse(impulse_idx, 
                                         impulse_idx).topRows(dimv_)
          = kkt_matrix.impulse[impulse_idx].Pq().transpose();
      constraint_factorization.T_impulse(
          impulse_idx, impulse_idx).bottomRows(dimv_).setZero();
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_idx = i - N_ - N_impulse;
      riccati_factorizer.aux[impulse_idx].forwardRiccatiRecursionParallel(
          kkt_matrix.aux[impulse_idx], kkt_residual.aux[impulse_idx], 
          exist_state_constraint);
    }
    else {
      const int lift_idx = i - N_ - 2 * N_impulse;
      riccati_factorizer.lift[lift_idx].forwardRiccatiRecursionParallel(
          kkt_matrix.lift[lift_idx], kkt_residual.lift[lift_idx], 
          exist_state_constraint);
    }
  }
}


void RiccatiRecursion::forwardStateConstraintFactorization(
    RiccatiFactorizer& riccati_factorizer, 
    const OCPDiscretizer& ocp_discretizer, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, 
    RiccatiFactorization& riccati_factorization) {
  const bool exist_state_constraint = ocp_discretizer.existImpulse();
  riccati_factorizer[0].forwardStateConstraintFactorizationInitial(
      riccati_factorization[0]);
  for (int i=0; i<N_; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_idx = ocp_discretizer.impulseIndex(i);
      riccati_factorizer[i].forwardStateConstraintFactorization(
          riccati_factorization[i], kkt_matrix[i], kkt_residual[i],
          ocp_discretizer.dtau(i), riccati_factorization.impulse[impulse_idx],
          exist_state_constraint);
      riccati_factorizer.impulse[impulse_idx].forwardStateConstraintFactorization(
          riccati_factorization.impulse[impulse_idx],
          kkt_matrix.impulse[impulse_idx], kkt_residual.impulse[impulse_idx],
          riccati_factorization.aux[impulse_idx], exist_state_constraint);
      riccati_factorizer.aux[impulse_idx].forwardStateConstraintFactorization(
          riccati_factorization.aux[impulse_idx],
          kkt_matrix.aux[impulse_idx], kkt_residual.aux[impulse_idx],
          ocp_discretizer.dtau_aux(impulse_idx), 
          riccati_factorization[i+1], exist_state_constraint);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_idx = ocp_discretizer.liftIndex(i);
      riccati_factorizer[i].forwardStateConstraintFactorization(
          riccati_factorization[i], kkt_matrix[i], kkt_residual[i], 
          ocp_discretizer.dtau(i), riccati_factorization.lift[lift_idx], 
          exist_state_constraint);
      riccati_factorizer.lift[lift_idx].forwardStateConstraintFactorization(
          riccati_factorization.lift[lift_idx], kkt_matrix.lift[lift_idx], 
          kkt_residual.lift[lift_idx], ocp_discretizer.dtau_lift(lift_idx), 
          riccati_factorization[i+1], exist_state_constraint);
    }
    else {
      riccati_factorizer[i].forwardStateConstraintFactorization(
        riccati_factorization[i], kkt_matrix[i], kkt_residual[i], 
        ocp_discretizer.dtau(i), riccati_factorization[i+1], 
        exist_state_constraint);
    }
  }
}


void RiccatiRecursion::backwardStateConstraintFactorization(
    const RiccatiFactorizer& riccati_factorizer, 
    const OCPDiscretizer& ocp_discretizer, const KKTMatrix& kkt_matrix, 
    StateConstraintRiccatiFactorization& constraint_factorization) const {
  const int num_constraint = ocp_discretizer.numImpulseStages();
  #pragma omp parallel for num_threads(nproc_)
  for (int constraint_idx=0; constraint_idx<num_constraint; ++constraint_idx) {
    const int time_stage_before_constraint 
        = ocp_discretizer.timeStageBeforeImpulse(constraint_idx);
    riccati_factorizer[time_stage_before_constraint].backwardStateConstraintFactorization(
        constraint_factorization.T_impulse(constraint_idx, constraint_idx), 
        kkt_matrix[time_stage_before_constraint], 
        ocp_discretizer.dtau(time_stage_before_constraint), 
        constraint_factorization.T(constraint_idx, time_stage_before_constraint));
    for (int i=time_stage_before_constraint-1; i>=0; --i) {
      if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
        const int impulse_idx = ocp_discretizer.impulseIndex(i);
        riccati_factorizer.aux[impulse_idx].backwardStateConstraintFactorization(
            constraint_factorization.T(constraint_idx, i+1), 
            kkt_matrix.aux[impulse_idx], ocp_discretizer.dtau_aux(impulse_idx), 
            constraint_factorization.T_aux(constraint_idx, impulse_idx));
        riccati_factorizer.impulse[impulse_idx].backwardStateConstraintFactorization(
            constraint_factorization.T_aux(constraint_idx, impulse_idx), 
            kkt_matrix.impulse[impulse_idx], 
            constraint_factorization.T_impulse(constraint_idx, impulse_idx));
        riccati_factorizer[i].backwardStateConstraintFactorization(
            constraint_factorization.T_impulse(constraint_idx, impulse_idx), 
            kkt_matrix[i], ocp_discretizer.dtau(i), 
            constraint_factorization.T(constraint_idx, i));
      }
      else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
        const int lift_idx = ocp_discretizer.liftIndex(i);
        riccati_factorizer.lift[lift_idx].backwardStateConstraintFactorization(
            constraint_factorization.T(constraint_idx, i+1), 
            kkt_matrix.lift[lift_idx], ocp_discretizer.dtau_lift(lift_idx), 
            constraint_factorization.T_lift(constraint_idx, lift_idx));
        riccati_factorizer[i].backwardStateConstraintFactorization(
            constraint_factorization.T_lift(constraint_idx, lift_idx), 
            kkt_matrix[i], ocp_discretizer.dtau(i), 
            constraint_factorization.T(constraint_idx, i));
      }
      else {
        riccati_factorizer[i].backwardStateConstraintFactorization(
            constraint_factorization.T(constraint_idx, i+1), 
            kkt_matrix[i], ocp_discretizer.dtau(i), 
            constraint_factorization.T(constraint_idx, i));
      }
    }
  }
}


void RiccatiRecursion::forwardRiccatiRecursion(
    const RiccatiFactorizer& riccati_factorizer, 
    const OCPDiscretizer& ocp_discretizer, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, 
    const RiccatiFactorization& riccati_factorization, Direction& d) {
  const bool exist_state_constraint = ocp_discretizer.existImpulse();
  for (int i=0; i<N_; ++i) {
    if (ocp_discretizer.isTimeStageBeforeImpulse(i)) {
      const int impulse_idx = ocp_discretizer.impulseIndex(i);
      riccati_factorizer[i].forwardRiccatiRecursion(
          kkt_matrix[i], kkt_residual[i], 
          riccati_factorization.impulse[impulse_idx], d[i], 
          ocp_discretizer.dtau(i), d.impulse[impulse_idx], 
          exist_state_constraint);
      riccati_factorizer.impulse[impulse_idx].forwardRiccatiRecursion(
          kkt_matrix.impulse[impulse_idx], kkt_residual.impulse[impulse_idx],
          d.impulse[impulse_idx], d.aux[impulse_idx]);
      riccati_factorizer.aux[impulse_idx].forwardRiccatiRecursion(
          kkt_matrix.aux[impulse_idx], kkt_residual.aux[impulse_idx],
          riccati_factorization[i+1], d.aux[impulse_idx], 
          ocp_discretizer.dtau_aux(impulse_idx), d[i+1], exist_state_constraint);
    }
    else if (ocp_discretizer.isTimeStageBeforeLift(i)) {
      const int lift_idx = ocp_discretizer.liftIndex(i);
      riccati_factorizer[i].forwardRiccatiRecursion(
          kkt_matrix[i], kkt_residual[i], riccati_factorization.lift[lift_idx], 
          d[i], ocp_discretizer.dtau(i), d.lift[lift_idx], 
          exist_state_constraint);
      riccati_factorizer.lift[lift_idx].forwardRiccatiRecursion(
          kkt_matrix.lift[lift_idx], kkt_residual.lift[lift_idx],
          riccati_factorization[i+1], d.lift[lift_idx], 
          ocp_discretizer.dtau_lift(lift_idx), d[i+1], exist_state_constraint);
    }
    else {
      riccati_factorizer[i].forwardRiccatiRecursion(
          kkt_matrix[i], kkt_residual[i], riccati_factorization[i+1], 
          d[i], ocp_discretizer.dtau(i), d[i+1], exist_state_constraint);
    }
  }
}

} // namespace idocp