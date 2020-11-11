#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HXX_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"

#include <omp.h>
#include <cassert>

namespace idocp {

inline StateConstraintRiccatiFactorizer::StateConstraintRiccatiFactorizer(
    const Robot& robot, const int max_num_impulse, const int nproc) 
  : llts_(max_num_impulse, Eigen::LLT<Eigen::MatrixXd>()),
    max_num_impulse_(max_num_impulse),
    dimv_(robot.dimv()),
    nproc_(nproc) { 
}


inline StateConstraintRiccatiFactorizer::StateConstraintRiccatiFactorizer() 
  : llts_(),
    max_num_impulse_(0),
    dimv_(0),
    nproc_(0) { 
}


inline StateConstraintRiccatiFactorizer::~StateConstraintRiccatiFactorizer() { 
}


template <typename VectorType>
inline void StateConstraintRiccatiFactorizer::computeLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence,
    const std::vector<RiccatiFactorization>& impulse_riccati_factorization,
    std::vector<StateConstraintRiccatiFactorization>& constraint_factorization,
    const Eigen::MatrixBase<VectorType>& dx0,
    std::vector<ImpulseSplitDirection>& d) {
  assert(impulse_riccati_factorization.size() == max_num_impulse_);
  assert(constraint_factorization.size() == max_num_impulse_);
  assert(d.size() == max_num_impulse_);
  const int num_impulse = contact_sequence.numImpulse();
  assert(num_impulse <= max_num_impulse_);
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<num_impulse; ++i) {
    factorizeLinearProblem(impulse_riccati_factorization[i], 
                           constraint_factorization[i], dx0);
    llts_[i].compute(constraint_factorization[i].ENEt());
    assert(llts_[i].info() == Eigen::Success);
  }
  for (int i=num_impulse-1; i>=0; --i) {
    for (int j=i+1; j<num_impulse; ++j) {
      constraint_factorization[i].e().noalias() 
          -= constraint_factorization[i].EN() 
              * constraint_factorization[i].T_impulse(j) * d[j].dxi();
    }
    d[i].dxi() = llts_[i].solve(constraint_factorization[i].e());
  }
}


template <typename VectorType>
inline void StateConstraintRiccatiFactorizer::factorizeLinearProblem(
    const RiccatiFactorization& impulse_riccati_factorization,
    StateConstraintRiccatiFactorization& constraint_factorization,
    const Eigen::MatrixBase<VectorType>& dx0) {
  constraint_factorization.EN().noalias() 
      = constraint_factorization.Eq() 
          * impulse_riccati_factorization.N.topRows(dimv_);
  constraint_factorization.ENEt().noalias() 
      = constraint_factorization.ENq() 
          * constraint_factorization.Eq().transpose();
  constraint_factorization.e().noalias() 
      += constraint_factorization.Eq() 
          * (impulse_riccati_factorization.Pi.topRows(dimv_) * dx0 
              + impulse_riccati_factorization.pi.head(dimv_)); 
}

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZER_HXX_ 