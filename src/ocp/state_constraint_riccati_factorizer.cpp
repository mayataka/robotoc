#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"

#include <omp.h>
#include <cassert>

namespace idocp {

StateConstraintRiccatiFactorizer::StateConstraintRiccatiFactorizer(
    const Robot& robot, const int N, const int max_num_impulse, const int nproc) 
  : ldlt_(Eigen::LDLT<Eigen::MatrixXd>()),
    lp_factorizer_(max_num_impulse, StateConstraintRiccatiLPFactorizer(robot)),
    N_(N),
    max_num_impulse_(max_num_impulse),
    nproc_(nproc) { 
}


StateConstraintRiccatiFactorizer::StateConstraintRiccatiFactorizer() 
  : ldlt_(),
    lp_factorizer_(),
    N_(0),
    max_num_impulse_(0),
    nproc_(0) { 
}


StateConstraintRiccatiFactorizer::~StateConstraintRiccatiFactorizer() { 
}


void StateConstraintRiccatiFactorizer::computeLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence,
    const HybridRiccatiFactorization& riccati_factorization,
    StateConstraintRiccatiFactorization& constraint_factorization,
    Direction& d) {
  assert(riccati_factorization.impulse.size() == max_num_impulse_);
  assert(d.impulse.size() == max_num_impulse_);
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  assert(num_impulse <= max_num_impulse_);
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<num_impulse; ++i) {
    lp_factorizer_[i].factorizeLinearProblem(contact_sequence, 
                                             riccati_factorization.impulse[i], 
                                             constraint_factorization, 
                                             d[0].dx(), i);
  }
  constraint_factorization.ENT().triangularView<Eigen::StrictlyLower>() 
      = constraint_factorization.ENT().transpose().triangularView<Eigen::StrictlyLower>();
  ldlt_.compute(constraint_factorization.ENT());
  assert(ldlt_.info() == Eigen::Success);
  constraint_factorization.dxi() = ldlt_.solve(constraint_factorization.e());
  for (int i=0; i<num_impulse; ++i) {
    d.impulse[i].dxi() = constraint_factorization.dxi(i);
  }
}


void StateConstraintRiccatiFactorizer::aggregateLagrangeMultiplierDirection(
    const StateConstraintRiccatiFactorization& constraint_factorization,
    const ContactSequence& contact_sequence, const Direction& d,
    HybridRiccatiFactorization& riccati_factorization) const {
  const int N_impulse = contact_sequence.totalNumImpulseStages();
  const int N_lift = contact_sequence.totalNumLiftStages();
  const int N_all = N_ + 2*N_impulse + N_lift;
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<N_all; ++i) {
    if (i < N_) {
      riccati_factorization[i].n.setZero();
      for (int constraint_idx=N_impulse-1; constraint_idx>=0; --constraint_idx) {
        if (contact_sequence.timeStageBeforeImpulse(constraint_idx) < i) {
          break;
        }
        else {
          riccati_factorization[i].n.noalias()
              += constraint_factorization.T(constraint_idx, i) 
                  * d.impulse[constraint_idx].dxi();
        }
      }
    } 
    else if (i < N_+N_impulse) {
      const int impulse_idx = i - N_;
      riccati_factorization.impulse[impulse_idx].n.setZero();
      for (int constraint_idx=N_impulse-1; constraint_idx>=impulse_idx; --constraint_idx) {
        riccati_factorization.impulse[impulse_idx].n.noalias()
            += constraint_factorization.T_impulse(constraint_idx, impulse_idx) 
                * d.impulse[constraint_idx].dxi();
      }
    }
    else if (i < N_+2*N_impulse) {
      const int impulse_idx = i - N_ - N_impulse;
      riccati_factorization.aux[impulse_idx].n.setZero();
      for (int constraint_idx=N_impulse-1; constraint_idx>impulse_idx; --constraint_idx) {
        riccati_factorization.aux[impulse_idx].n.noalias()
            += constraint_factorization.T_aux(constraint_idx, impulse_idx) 
                * d.impulse[constraint_idx].dxi();
      }
    }
    else {
      const int lift_idx = i - N_ - 2 * N_impulse;
      const int time_stage_before_lift 
          = contact_sequence.timeStageBeforeLift(lift_idx);
      riccati_factorization.lift[lift_idx].n.setZero();
      for (int constraint_idx=N_impulse-1; constraint_idx>=0; --constraint_idx) {
        if (contact_sequence.timeStageBeforeImpulse(constraint_idx) < time_stage_before_lift) {
          break;
        }
        else {
          riccati_factorization.lift[lift_idx].n.noalias()
              += constraint_factorization.T_lift(constraint_idx, lift_idx) 
                  * d.impulse[constraint_idx].dxi();
        }
      }
    }
  }
}

} // namespace idocp