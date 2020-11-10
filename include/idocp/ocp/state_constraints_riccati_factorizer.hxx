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
    dimf_(0),
    nproc_(nproc),
    is_active_(false) { 
}


inline StateConstraintsRiccatiFactorizer::StateConstraintsRiccatiFactorizer() 
  : llts_(),
    N_(0),
    dimv_(0),
    dimx_(0),
    dimf_(0), 
    nproc_(0),
    is_active_(false) { 
}


inline StateConstraintsRiccatiFactorizer::~StateConstraintsRiccatiFactorizer() { 
}


inline void StateConstraintsRiccatiFactorizer::computeLagrangeMultiplierDirection(
    const ContactSequence& contact_sequence,
    const std::vector<StateConstraintRiccatiFactorization>& factorizations,
    const std::vector<RiccatiSolution>& riccati_solutions,
    std::vector<ImpulseSplitDirection>& d) {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<N_; ++i) {
    if (contact_sequence_.existImpulse(i)) {
      llts_[i].compute(factorizations[i].ENEt());
    }
  }
}

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINTS_RICCATI_FACTORIZER_HXX_ 