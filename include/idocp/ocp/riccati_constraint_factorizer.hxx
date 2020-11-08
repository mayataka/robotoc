#ifndef IDOCP_RICCATI_CONSTRAINT_FACTORIZER_HXX_
#define IDOCP_RICCATI_CONSTRAINT_FACTORIZER_HXX_

#include "idocp/ocp/riccati_constraint_factorizer.hpp"

#include <omp.h>


namespace idocp {

inline RiccatiConstraintFactorizer::RiccatiConstraintFactorizer(
    const Robot& robot, const int N, const int nproc)
  : llts_(max_constrained_steps, Eigen::LLT<Eigen::MatrixXd>()),
    N_(N),
    nproc_(nproc) {
}


inline RiccatiConstraintFactorizer::RiccatiConstraintFactorizer()
  : llts_(),
    N_(0),
    nproc_(0) {
}


inline void RiccatiConstraintFactorizer::solve(
    const std::vector<RiccatiSolution>& riccatis) {

}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
inline void RiccatiConstraintFactorizer::computeENEtinv(
    Eigen::LLT<Eigen::MatrixXd>& llt,
    const Eigen::MatrixBase<MatrixType1>& Eq, 
    const Eigen::MatrixBase<MatrixType2>& Nqq, 
    const Eigen::MatrixBase<MatrixType3>& ENEtinv) {
  assert(Eq.rows() == dimf_);
  assert(Eq.cols() == dimv_);
  assert(Nqq.rows() == dimv_);
  assert(Nqq.cols() == dimv_);
  assert(ENEtinv.rows() == dimf_);
  assert(ENEtinv.cols() == dimf_);
  EqNqq_.topRows(dimf_).noalias() = Eq * Nqq;
  ENEt_.topLeftCorner(dimf_, dimf_).noalias() 
      = EqNqq_.topRows(dimf_) * Eq.transpose();
  llt.compute(ENEt_.topLeftCorner(dimf_, dimf_));
  assert(llt.info() == Eigen::Success);
  const_cast<Eigen::MatrixBase<MatrixType3>&>(ENEtinv) 
      = llt.solve(Eigen::MatrixXd::Identity(dimf_, dimf_));
}


inline void RiccatiConstraintFactorizer::invertBlockLowerTriangular(
    const std::vector<RiccatiSolution>& riccati_impulse, 
    const std::vector<RiccatiConstraintSolution>& riccati_constraint_impulse) {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<riccatis.size(); ++i) {
    if (riccatis.isActive()) {
      riccatis_constraint_impulse[i].computeENEt(
          riccati_constraint_impulseT().topRows(dimv_), 
                              riccatis[i].N.topLeftCorner(dimv_, dimv_));
      llts_[i].compute()
      riccatis.computeENEt();
      riccatis.solve();
      riccatis.llt.compute(riccatis[i].)
    }
  }

  invertBlockDiagonal();
}


inline void RiccatiConstraintFactorizer::invertBlockDiagonal() {
  #pragma omp parallel for num_threads(num_proc_)
  for (int i=0; i<active_constriants; ++i) {
      computeENEtinv();
  }
}

} // namespace idocp

#endif // IDOCP_RICCATI_CONSTRAINT_FACTORIZER_HXX_