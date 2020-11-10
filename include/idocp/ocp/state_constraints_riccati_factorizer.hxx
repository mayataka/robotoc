#ifndef IDOCP_STATE_CONSTRAINTS_RICCATI_FACTORIZER_HXX_
#define IDOCP_STATE_CONSTRAINTS_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/state_constraints_riccati_factorizer.hpp"

#include <omp.h>
#include <cassert>

namespace idocp {

inline StateConstraintsRiccatiFactorizer::StateConstraintsRiccatiFactorizer(
    const Robot& robot, const int N, const int nproc) 
  : llts_(N, Eigen::LLT<Eigen::MatrixXd>()),
    N_(N),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    nproc_(nproc) { 
}


inline StateConstraintsRiccatiFactorizer::StateConstraintsRiccatiFactorizer() 
  : llts_(),
    N_(0),
    dimv_(0),
    dimx_(0),
    nproc_(0) { 
}


inline StateConstraintsRiccatiFactorizer::~StateConstraintsRiccatiFactorizer() { 
}


template <typename VectorType>
inline void StateConstraintsRiccatiFactorizer::computeLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence,
    std::vector<StateConstraintRiccatiFactorization>& factorizations,
    const std::vector<RiccatiSolution>& riccati_solutions,
    const Eigen::MatrixBase<VectorType>& dx0,
    std::vector<ImpulseSplitDirection>& d) {
  assert(factorizations.size() == N_);
  assert(riccati_solutions.size() == N_);
  assert(d.size() == N_);
  #pragma omp parallel for num_threads(nproc_)
  for (int i=0; i<N_; ++i) {
    if (contact_sequence.existImpulse(i)) {
      llts_[i].compute(factorizations[i].ENEt());
      factorizations[i].e().noalias() += factorizations[i].Eq() * (riccati_solutions[i].Pi * dx0 + riccati_solutions[i].pi).topRows(dimv_); 
    }
  }
  for (int i=N_-1; i>=0; --i) {
    for (int j=i; j<N_; ++j) {
      if (contact_sequence.existImpulse(i)) {
        // factorizations.rhs().noalias() -= E * N * T * d[j].xi_stack();
      }
    }
    d[i].dxi() = llts_[i].solve(factorizations[i].e());
  }
}

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINTS_RICCATI_FACTORIZER_HXX_ 