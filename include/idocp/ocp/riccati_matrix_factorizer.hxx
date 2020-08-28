#ifndef IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_
#define IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_

#include "Eigen/LU"

#include <assert.h>

namespace idocp {

inline RiccatiMatrixFactorizer::RiccatiMatrixFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    dsubtract_dq_(),
    dsubtract_dq_prev_inv_() {
  if (robot.has_floating_base()) {
    dsubtract_dq_.resize(kDimFloatingBase, kDimFloatingBase);
    dsubtract_dq_.setZero();
    dsubtract_dq_prev_inv_.resize(kDimFloatingBase, kDimFloatingBase);
    dsubtract_dq_prev_inv_.setZero();
  }
}


inline RiccatiMatrixFactorizer::RiccatiMatrixFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dsubtract_dq_(),
    dsubtract_dq_prev_inv_() {
}


inline RiccatiMatrixFactorizer::~RiccatiMatrixFactorizer() {
}


template <typename MatrixType>
inline void RiccatiMatrixFactorizer::setStateEquationDerivative(
    const Eigen::MatrixBase<MatrixType>& dsubtract_dq) {
  assert(dsubtract_dq.rows() >= kDimFloatingBase);
  assert(dsubtract_dq.cols() >= kDimFloatingBase);
  if (has_floating_base_) {
    dsubtract_dq_ 
        = dsubtract_dq.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>();
  }
}


template <typename MatrixType>
inline void RiccatiMatrixFactorizer::setStateEquationDerivativeInverse(
    const Eigen::MatrixBase<MatrixType>& dsubtract_dq_prev) {
  assert(dsubtract_dq_prev.rows() >= kDimFloatingBase);
  assert(dsubtract_dq_prev.cols() >= kDimFloatingBase);
  if (has_floating_base_) {
    dsubtract_dq_prev_inv_ 
            = - dsubtract_dq_prev.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>().inverse();
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4, typename MatrixType5, typename MatrixType6, 
          typename MatrixType7, typename MatrixType8>
inline void RiccatiMatrixFactorizer::factorizeF(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Pqq_next,
    const Eigen::MatrixBase<MatrixType2>& Pqv_next,
    const Eigen::MatrixBase<MatrixType3>& Pvq_next,
    const Eigen::MatrixBase<MatrixType4>& Pvv_next,
    const Eigen::MatrixBase<MatrixType5>& Qqq,
    const Eigen::MatrixBase<MatrixType6>& Qqv,
    const Eigen::MatrixBase<MatrixType7>& Qvq,
    const Eigen::MatrixBase<MatrixType8>& Qvv) const {
  assert(dtau > 0);
  assert(Pqq_next.rows() == dimv_);
  assert(Pqq_next.cols() == dimv_);
  assert(Pqv_next.rows() == dimv_);
  assert(Pqv_next.cols() == dimv_);
  assert(Pvq_next.rows() == dimv_);
  assert(Pvq_next.cols() == dimv_);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Qqq.rows() == dimv_);
  assert(Qqq.cols() == dimv_);
  assert(Qqv.rows() == dimv_);
  assert(Qqv.cols() == dimv_);
  assert(Qvq.rows() == dimv_);
  assert(Qvq.cols() == dimv_);
  assert(Qvv.rows() == dimv_);
  assert(Qvv.cols() == dimv_);
  if (has_floating_base_) {
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).template topLeftCorner<kDimFloatingBase, kDimFloatingBase>().noalias() 
        += dsubtract_dq_.transpose() 
            * Pqq_next.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>()
            * dsubtract_dq_;
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).topRightCorner(kDimFloatingBase, dimv_-kDimFloatingBase).noalias()
        += dsubtract_dq_.transpose()
            * Pqq_next.topRightCorner(kDimFloatingBase, dimv_-kDimFloatingBase);
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).bottomLeftCorner(dimv_-kDimFloatingBase, kDimFloatingBase).noalias()
        += Pqq_next.bottomLeftCorner(dimv_-kDimFloatingBase, kDimFloatingBase) 
            * dsubtract_dq_;
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).bottomRightCorner(dimv_-kDimFloatingBase, dimv_-kDimFloatingBase).noalias()
        += Pqq_next.bottomRightCorner(dimv_-kDimFloatingBase, dimv_-kDimFloatingBase);

    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).template topRows<kDimFloatingBase>().noalias() 
        += dtau * dsubtract_dq_.transpose() 
                * Pqq_next.template topRows<kDimFloatingBase>();
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).bottomRows(dimv_-kDimFloatingBase).noalias() 
        += dtau * Pqq_next.bottomRows(dimv_-kDimFloatingBase);
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).template topRows<kDimFloatingBase>().noalias() 
        += dsubtract_dq_.transpose() 
                * Pqv_next.template topRows<kDimFloatingBase>();
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).bottomRows(dimv_-kDimFloatingBase).noalias() 
        += Pqv_next.bottomRows(dimv_-kDimFloatingBase);

    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).template leftCols<kDimFloatingBase>().noalias() 
        += dtau * Pqq_next.template leftCols<kDimFloatingBase>()
                * dsubtract_dq_;
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).rightCols(dimv_-kDimFloatingBase).noalias() 
        += dtau * Pqq_next.rightCols(dimv_-kDimFloatingBase);
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).template leftCols<kDimFloatingBase>().noalias() 
        += Pvq_next.template leftCols<kDimFloatingBase>() * dsubtract_dq_;
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).rightCols(dimv_-kDimFloatingBase).noalias() 
        += Pvq_next.rightCols(dimv_-kDimFloatingBase);
  }
  else {
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).noalias() += Pqq_next;
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).noalias() += dtau * Pqq_next;
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).noalias() += Pqv_next;
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).noalias() += dtau * Pqq_next;
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).noalias() += Pvq_next;
  }
  (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).noalias() += (dtau*dtau) * Pqq_next;
  (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).noalias() += dtau * (Pqv_next + Pvq_next);
  (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).noalias() += Pvv_next;
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4>
inline void RiccatiMatrixFactorizer::factorizeH(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Pqv_next, 
    const Eigen::MatrixBase<MatrixType2>& Pvv_next, 
    const Eigen::MatrixBase<MatrixType3>& Qqa, 
    const Eigen::MatrixBase<MatrixType4>& Qva) const {
  assert(dtau > 0);
  assert(Pqv_next.rows() == dimv_);
  assert(Pqv_next.cols() == dimv_);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Qqa.rows() == dimv_);
  assert(Qqa.cols() == dimv_);
  assert(Qva.rows() == dimv_);
  assert(Qva.cols() == dimv_);
  if (has_floating_base_) {
    (const_cast<Eigen::MatrixBase<MatrixType3>&>(Qqa)).template topRows<kDimFloatingBase>().noalias()
        += dtau * dsubtract_dq_.transpose()
                * Pqv_next.template topRows<kDimFloatingBase>();
    (const_cast<Eigen::MatrixBase<MatrixType3>&>(Qqa)).bottomRows(dimv_-kDimFloatingBase).noalias()
        += dtau * Pqv_next.bottomRows(dimv_-kDimFloatingBase);
  }
  else {
    (const_cast<Eigen::MatrixBase<MatrixType3>&>(Qqa)).noalias() += dtau * Pqv_next;
  }
  (const_cast<Eigen::MatrixBase<MatrixType4>&>(Qva)).noalias() += (dtau*dtau) * Pqv_next;
  (const_cast<Eigen::MatrixBase<MatrixType4>&>(Qva)).noalias() += dtau * Pvv_next;
}


template <typename MatrixType1, typename MatrixType2>
inline void RiccatiMatrixFactorizer::factorizeG(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Pvv_next, 
    const Eigen::MatrixBase<MatrixType2>& Qaa) const {
  assert(dtau > 0);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Qaa.rows() == dimv_);
  assert(Qaa.cols() == dimv_);
  (const_cast<Eigen::MatrixBase<MatrixType2>&>(Qaa)).noalias()
      += (dtau*dtau) * Pvv_next;
}


template <typename MatrixType1, typename MatrixType2, typename VectorType1,
          typename VectorType2, typename VectorType3, typename VectorType4>
inline void RiccatiMatrixFactorizer::factorize_la(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Pvq_next, 
    const Eigen::MatrixBase<MatrixType2>& Pvv_next, 
    const Eigen::MatrixBase<VectorType1>& Fq, 
    const Eigen::MatrixBase<VectorType2>& Fv, 
    const Eigen::MatrixBase<VectorType3>& sv_next, 
    const Eigen::MatrixBase<VectorType4>& la) const {
  assert(dtau > 0);
  assert(Pvq_next.rows() == dimv_);
  assert(Pvq_next.cols() == dimv_);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Fq.size() == dimv_);
  assert(Fv.size() == dimv_);
  assert(sv_next.size() == dimv_);
  assert(la.size() == dimv_);
  (const_cast<Eigen::MatrixBase<VectorType4>&>(la)).noalias() += dtau * Pvq_next * Fq; 
  (const_cast<Eigen::MatrixBase<VectorType4>&>(la)).noalias() += dtau * Pvv_next * Fv; 
  (const_cast<Eigen::MatrixBase<VectorType4>&>(la)).noalias() -= dtau * sv_next; 
}


template <typename MatrixType1, typename MatrixType2>
inline void RiccatiMatrixFactorizer::correctP(
    const Eigen::MatrixBase<MatrixType1>& Pqq, 
    const Eigen::MatrixBase<MatrixType2>& Pqv) const {
  assert(Pqq.rows() == dimv_);
  assert(Pqq.cols() == dimv_);
  assert(Pqv.rows() == dimv_);
  assert(Pqv.cols() == dimv_);
  if (has_floating_base_) {
    (const_cast<Eigen::MatrixBase<MatrixType1>&>(Pqq)).template topRows<kDimFloatingBase>()
        = dsubtract_dq_prev_inv_.transpose() * Pqq.template topRows<kDimFloatingBase>();
    (const_cast<Eigen::MatrixBase<MatrixType2>&>(Pqv)).template topRows<kDimFloatingBase>()
        = dsubtract_dq_prev_inv_.transpose() * Pqv.template topRows<kDimFloatingBase>();
  }
}


template <typename VectorType>
inline void RiccatiMatrixFactorizer::correct_s(
    const Eigen::MatrixBase<VectorType>& sq) const {
  assert(sq.size() == dimv_);
  if (has_floating_base_) {
    (const_cast<Eigen::MatrixBase<VectorType>&>(sq)).template head<kDimFloatingBase>()
        = dsubtract_dq_prev_inv_.transpose() * sq.template head<kDimFloatingBase>();
  }
}

} // namespace idocp

#endif // IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_