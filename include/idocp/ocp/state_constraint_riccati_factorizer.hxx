#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HXX_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"

#include <omp.h>
#include <cassert>

namespace idocp {

inline StateConstraintRiccatiFactorizer::StateConstraintRiccatiFactorizer(
    const Robot& robot, const int max_num_impulse, const int nproc) 
  : ldlt_(Eigen::LDLT<Eigen::MatrixXd>()),
    lp_factorizer_(max_num_impulse, StateConstraintRiccatiLPFactorizer(robot)),
    max_num_impulse_(max_num_impulse),
    nproc_(nproc) { 
}


inline StateConstraintRiccatiFactorizer::StateConstraintRiccatiFactorizer() 
  : ldlt_(),
    lp_factorizer_(),
    max_num_impulse_(0),
    nproc_(0) { 
}


inline StateConstraintRiccatiFactorizer::~StateConstraintRiccatiFactorizer() { 
}


template <typename VectorType>
inline void StateConstraintRiccatiFactorizer::computeLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence,
    const std::vector<RiccatiFactorization>& impulse_riccati_factorization,
    StateConstraintRiccatiFactorization& constraint_factorization,
    const Eigen::MatrixBase<VectorType>& dx0,
    std::vector<ImpulseSplitDirection>& d_impulse) {
  assert(impulse_riccati_factorization.size() == max_num_impulse_);
  assert(d_impulse.size() == max_num_impulse_);
  const int num_impulse = contact_sequence.totalNumImpulseStages();
  assert(num_impulse <= max_num_impulse_);
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<num_impulse; ++i) {
    lp_factorizer_[i].factorizeLinearProblem(contact_sequence, 
                                             impulse_riccati_factorization[i], 
                                             constraint_factorization, dx0, i);
  }
  constraint_factorization.ENT().triangularView<Eigen::StrictlyLower>() 
      = constraint_factorization.ENT().transpose().triangularView<Eigen::StrictlyLower>();
  ldlt_.compute(constraint_factorization.ENT());
  assert(ldlt_.info() == Eigen::Success);
  constraint_factorization.dxi() = ldlt_.solve(constraint_factorization.e());
  for (int i=0; i<num_impulse; ++i) {
    d_impulse[i].dxi() = constraint_factorization.dxi(i);
  }
}

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HXX_ 