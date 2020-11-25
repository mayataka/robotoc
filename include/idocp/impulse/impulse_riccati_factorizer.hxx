#ifndef IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_
#define IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_

#include "idocp/impulse/impulse_riccati_factorizer.hpp"

namespace idocp {

inline ImpulseRiccatiFactorizer::ImpulseRiccatiFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    backward_recursion_(robot),
    forward_recursion_(robot) {
}


inline ImpulseRiccatiFactorizer::ImpulseRiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    backward_recursion_(),
    forward_recursion_() {
}


inline ImpulseRiccatiFactorizer::~ImpulseRiccatiFactorizer() {
}


inline void ImpulseRiccatiFactorizer::backwardRiccatiRecursion(
    const RiccatiFactorization& riccati_next, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual, RiccatiFactorization& riccati) {
  backward_recursion_.factorizeKKTMatrix(riccati_next, kkt_matrix);
  backward_recursion_.factorizeRiccatiFactorization(riccati_next, kkt_matrix, 
                                                    kkt_residual, riccati);
  // const int dimv = dimv_;
  // const int dimx = 2 * dimv;
  // Eigen::MatrixXd P_next = Eigen::MatrixXd::Zero(dimx, dimx);
  // P_next.topLeftCorner(dimv, dimv) = riccati_next.Pqq;
  // P_next.topRightCorner(dimv, dimv) = riccati_next.Pqv;
  // P_next.bottomLeftCorner(dimv, dimv) = riccati_next.Pvq;
  // P_next.bottomRightCorner(dimv, dimv) = riccati_next.Pvv;
  // // std::cout << "Fxx impulse backward" << std::endl;
  // // std::cout << kkt_matrix.Fxx() << std::endl;
  // // kkt_matrix.Qxx() += kkt_matrix.Fxx().transpose() * P_next * kkt_matrix.Fxx();
  // assert(kkt_matrix.Qxx().isApprox(kkt_matrix.Fxx().transpose() * P_next * kkt_matrix.Fxx()));
  // const Eigen::MatrixXd P = kkt_matrix.Qxx();
  // riccati.Pqq = P.topLeftCorner(dimv, dimv);
  // riccati.Pqv = P.topRightCorner(dimv, dimv);
  // riccati.Pvq = P.bottomLeftCorner(dimv, dimv);
  // riccati.Pvv = P.bottomRightCorner(dimv, dimv);
  // Eigen::VectorXd s_next = Eigen::VectorXd::Zero(dimx);
  // s_next.head(dimv) = riccati_next.sq;
  // s_next.tail(dimv) = riccati_next.sv;
  // const Eigen::VectorXd s = kkt_matrix.Fxx().transpose() * s_next
  //                           - kkt_matrix.Fxx().transpose() * P_next * kkt_residual.Fx()
  //                           - kkt_residual.lx();
  // riccati.sq = s.head(dimv);
  // riccati.sv = s.tail(dimv);
}


inline void ImpulseRiccatiFactorizer::forwardRiccatiRecursionSerial(
    const RiccatiFactorization& riccati, const ImpulseKKTMatrix& kkt_matrix, 
    const ImpulseKKTResidual& kkt_residual, 
    RiccatiFactorization& riccati_next) {
  forward_recursion_.factorizeStateTransition(riccati, kkt_matrix, kkt_residual, 
                                              riccati_next);
  forward_recursion_.factorizeStateConstraintFactorization(riccati, kkt_matrix,
                                                           riccati_next);
  // const int dimv = dimv_;
  // const int dimx = 2 * dimv;
  // std::cout << "Fxx impulse forward" << std::endl;
  // std::cout << kkt_matrix.Fxx() << std::endl;
  // std::cout << "Fx impulse forward" << std::endl;
  // std::cout << kkt_residual.Fx().transpose() << std::endl;
  // riccati_next.Pi = kkt_matrix.Fxx() * riccati.Pi;
  // riccati_next.pi = kkt_matrix.Fxx() * riccati.pi + kkt_residual.Fx();
}


template <typename MatrixType1, typename MatrixType2>
inline void ImpulseRiccatiFactorizer::backwardStateConstraintFactorization(
    const Eigen::MatrixBase<MatrixType1>& T_next,
    const ImpulseKKTMatrix& kkt_matrix, 
    const Eigen::MatrixBase<MatrixType2>& T) const {
  assert(T_next.rows() == T.rows());
  assert(T_next.rows() == T.rows());
  if (has_floating_base_) {
    const_cast<Eigen::MatrixBase<MatrixType2>&> (T).topRows(dimv_).noalias() 
        = kkt_matrix.Fqq().transpose() * T_next.topRows(dimv_);
  }
  else {
    const_cast<Eigen::MatrixBase<MatrixType2>&> (T).topRows(dimv_) 
        = T_next.topRows(dimv_);
  }
  const_cast<Eigen::MatrixBase<MatrixType2>&> (T).topRows(dimv_).noalias()
      += kkt_matrix.Fvq().transpose() * T_next.bottomRows(dimv_);
  const_cast<Eigen::MatrixBase<MatrixType2>&> (T).bottomRows(dimv_).noalias() 
      =  kkt_matrix.Fvv().transpose() * T_next.bottomRows(dimv_);
}


template <typename VectorType>
inline void ImpulseRiccatiFactorizer::computeStateDirection(
    const RiccatiFactorization& riccati, 
    const Eigen::MatrixBase<VectorType>& dx0, ImpulseSplitDirection& d) {
  d.dx().noalias() = riccati.Pi * dx0;
  d.dx().noalias() += riccati.pi;
  d.dx().noalias() -= riccati.N * riccati.n;
}


inline void ImpulseRiccatiFactorizer::computeCostateDirection(
    const RiccatiFactorization& riccati, ImpulseSplitDirection& d) {
  d.dlmd().noalias() = riccati.Pqq * d.dq();
  d.dlmd().noalias() += riccati.Pqv * d.dv();
  d.dlmd().noalias() -= riccati.sq;
  d.dgmm().noalias() = riccati.Pqv.transpose() * d.dq();
  d.dgmm().noalias() += riccati.Pvv * d.dv();
  d.dgmm().noalias() -= riccati.sv;
  d.dlmdgmm().noalias() += riccati.n;
}

} // namespace idocp

#endif // IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_