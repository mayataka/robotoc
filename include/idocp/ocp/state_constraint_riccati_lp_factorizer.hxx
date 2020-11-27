#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_LP_FACTORIZATER_HXX_ 
#define IDOCP_STATE_CONSTRAINT_RICCATI_LP_FACTORIZATER_HXX_ 

#include "idocp/ocp/state_constraint_riccati_lp_factorizer.hpp"

#include <cassert>

namespace idocp {

inline StateConstraintRiccatiLPFactorizer::StateConstraintRiccatiLPFactorizer(
    const Robot& robot) 
  : dimv_(robot.dimv()),
    Pidx0q_(Eigen::VectorXd::Zero(robot.dimv())) { 
}


inline StateConstraintRiccatiLPFactorizer::StateConstraintRiccatiLPFactorizer() 
  : dimv_(0),
    Pidx0q_() { 
}


inline StateConstraintRiccatiLPFactorizer::~StateConstraintRiccatiLPFactorizer() { 
}


template <typename VectorType>
inline void StateConstraintRiccatiLPFactorizer::factorizeLinearProblem(
    const RiccatiFactorization& impulse_riccati_factorization,
    StateConstraintRiccatiFactorization& constraint_factorization,
    const Eigen::MatrixBase<VectorType>& dx0, const int constraint_index) {
  constraint_factorization.EN(constraint_index).noalias() 
      = constraint_factorization.Eq(constraint_index) 
          * impulse_riccati_factorization.N.topRows(dimv_);
  constraint_factorization.ENEt(constraint_index).noalias() 
      = constraint_factorization.ENq(constraint_index) 
          * constraint_factorization.Eq(constraint_index).transpose();
  Pidx0q_.noalias() = impulse_riccati_factorization.Pi.topRows(dimv_) * dx0;
  Pidx0q_.noalias() += impulse_riccati_factorization.pi.head(dimv_);
  constraint_factorization.e(constraint_index).noalias() 
      += constraint_factorization.Eq(constraint_index) * Pidx0q_;
}

} // namespace idocp

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_LP_FACTORIZATER_HXX_ 