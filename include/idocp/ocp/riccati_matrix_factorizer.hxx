#ifndef IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_
#define IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_

namespace idocp {

inline RiccatiMatrixFactorizer::RiccatiMatrixFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    dintegrate_dq_(),
    dintegrate_dv_() {
  if (robot.has_floating_base()) {
    dintegrate_dq_.resize(robot.dimv(), robot.dimv());
    dintegrate_dv_.resize(robot.dimv(), robot.dimv());
  }
}


inline RiccatiMatrixFactorizer::RiccatiMatrixFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dintegrate_dq_(),
    dintegrate_dv_() {
}


inline RiccatiMatrixFactorizer::~RiccatiMatrixFactorizer() {
}


template <typename ConfigVectorType, typename TangentVectorType>
inline void RiccatiMatrixFactorizer::setIntegrationSensitivities(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<ConfigVectorType>& q,
    const Eigen::MatrixBase<TangentVectorType>& v) {
  assert(dtau > 0);
  assert(q.size() == robot.dimq());
  assert(v.size() == robot.dimv());
  if (has_floating_base_) {
    robot.dIntegrateConfiguration(q, v, dtau, dintegrate_dq_, dintegrate_dv_);
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
    const Eigen::MatrixBase<MatrixType8>& Qvv) {
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
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).template topLeftCorner<6, 6>().noalias() 
        += dintegrate_dq_.template topLeftCorner<6, 6>().transpose() 
            * Pqq_next.template topLeftCorner<6, 6>()
            * dintegrate_dq_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).topRightCorner(6, dimv_-6).noalias()
        += dintegrate_dq_.template topLeftCorner<6, 6>().transpose()
            * Pqq_next.topRightCorner(6, dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).bottomLeftCorner(dimv_-6, 6).noalias()
        += Pqq_next.bottomLeftCorner(dimv_-6, 6) 
            * dintegrate_dq_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);

    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).template topLeftCorner<6, 6>().noalias() 
        += dtau * dintegrate_dq_.template topLeftCorner<6, 6>().transpose() 
                * Pqq_next.template topLeftCorner<6, 6>() 
                * dintegrate_dv_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).topRightCorner(6, dimv_-6).noalias()
        += dtau * dintegrate_dq_.template topLeftCorner<6, 6>().transpose()
                * Pqq_next.topRightCorner(6, dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).bottomLeftCorner(dimv_-6, 6).noalias()
        += dtau * Pqq_next.bottomLeftCorner(dimv_-6, 6) 
                * dintegrate_dv_.topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += dtau * Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).template topRows<6>().noalias() 
        += dintegrate_dq_.template topLeftCorner<6, 6>().transpose()
            * Pqv_next.template topRows<6>();
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).bottomRows(dimv_-6).noalias() 
        += Pqv_next.bottomRows(dimv_-6);

    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).template topLeftCorner<6, 6>().noalias() 
        += dtau * dintegrate_dv_.template topLeftCorner<6, 6>().transpose() 
                * Pqq_next.template topLeftCorner<6, 6>() 
                * dintegrate_dq_.topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).topRightCorner(6, dimv_-6).noalias()
        += dtau *dintegrate_dv_.template topLeftCorner<6, 6>().transpose()
                * Pqq_next.topRightCorner(6, dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).bottomLeftCorner(dimv_-6, 6).noalias()
        += dtau * Pqq_next.bottomLeftCorner(dimv_-6, 6) 
                * dintegrate_dq_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += dtau * Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).template leftCols<6>().noalias() 
        += Pvq_next.template leftCols<6>() * dintegrate_dq_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).rightCols(dimv_-6).noalias() 
        += Pvq_next.rightCols(dimv_-6);

    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).template topLeftCorner<6, 6>().noalias() 
        += (dtau*dtau) * dintegrate_dv_.template topLeftCorner<6, 6>().transpose() 
                       * Pqq_next.template topLeftCorner<6, 6>()
                       * dintegrate_dv_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).topRightCorner(6, dimv_-6).noalias()
        += (dtau*dtau) * dintegrate_dv_.template topLeftCorner<6, 6>().transpose()
                       * Pqq_next.topRightCorner(6, dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).bottomLeftCorner(dimv_-6, 6).noalias()
        += (dtau*dtau) * Pqq_next.bottomLeftCorner(dimv_-6, 6) 
                       * dintegrate_dv_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).bottomRightCorner(dimv_-6, dimv_-6).noalias()
        += (dtau*dtau) * Pqq_next.bottomRightCorner(dimv_-6, dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).template leftCols<6>().noalias() 
        += dtau * Pvq_next.template leftCols<6>() * dintegrate_dv_.template topLeftCorner<6, 6>();
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).rightCols(dimv_-6).noalias() 
        += dtau * Pvq_next.rightCols(dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).template topRows<6>().noalias() 
        += dtau * dintegrate_dv_.template topLeftCorner<6, 6>().transpose()
                * Pqv_next.template topRows<6>();
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).bottomRows(dimv_-6).noalias() 
        += dtau * Pqv_next.bottomRows(dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).noalias() += Pvv_next;
  }
  else {
    (const_cast<Eigen::MatrixBase<MatrixType5>&>(Qqq)).noalias() += Pqq_next;
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).noalias() += dtau * Pqq_next;
    (const_cast<Eigen::MatrixBase<MatrixType6>&>(Qqv)).noalias() += Pqv_next;
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).noalias() += dtau * Pqq_next;
    (const_cast<Eigen::MatrixBase<MatrixType7>&>(Qvq)).noalias() += Pvq_next;
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).noalias() += (dtau*dtau) * Pqq_next;
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).noalias() += dtau * (Pqv_next + Pvq_next);
    (const_cast<Eigen::MatrixBase<MatrixType8>&>(Qvv)).noalias() += Pvv_next;
  }
}


template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
          typename MatrixType4>
inline void RiccatiMatrixFactorizer::factorizeH(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Pqv_next, 
    const Eigen::MatrixBase<MatrixType2>& Pvv_next, 
    const Eigen::MatrixBase<MatrixType3>& Qqa, 
    const Eigen::MatrixBase<MatrixType4>& Qva) {
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
    (const_cast<Eigen::MatrixBase<MatrixType3>&>(Qqa)).template topRows<6>().noalias()
        += dtau * dintegrate_dq_.template topLeftCorner<6, 6>().transpose()
                * Pqv_next.template topRows<6>();
    (const_cast<Eigen::MatrixBase<MatrixType3>&>(Qqa)).bottomRows(dimv_-6).noalias()
        += dtau * Pqv_next.bottomRows(dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType4>&>(Qva)).template topRows<6>().noalias()
        += (dtau*dtau) * dintegrate_dv_.template topLeftCorner<6, 6>().transpose()
                       * Pqv_next.template topRows<6>();
    (const_cast<Eigen::MatrixBase<MatrixType4>&>(Qva)).bottomRows(dimv_-6).noalias()
        += (dtau*dtau) * Pqv_next.bottomRows(dimv_-6);
    (const_cast<Eigen::MatrixBase<MatrixType4>&>(Qva)).noalias()
        += dtau * Pvv_next;
  }
  else {
    (const_cast<Eigen::MatrixBase<MatrixType3>&>(Qqa)).noalias() += dtau * Pqv_next;
    (const_cast<Eigen::MatrixBase<MatrixType4>&>(Qva)).noalias() += (dtau*dtau) * Pqv_next;
    (const_cast<Eigen::MatrixBase<MatrixType4>&>(Qva)).noalias() += dtau * Pvv_next;
  }
}


template <typename MatrixType1, typename MatrixType2>
inline void RiccatiMatrixFactorizer::factorizeG(
    const double dtau, const Eigen::MatrixBase<MatrixType1>& Pvv_next, 
    const Eigen::MatrixBase<MatrixType2>& Qaa) {
  assert(dtau > 0);
  assert(Pvv_next.rows() == dimv_);
  assert(Pvv_next.cols() == dimv_);
  assert(Qaa.rows() == dimv_);
  assert(Qaa.cols() == dimv_);
  (const_cast<Eigen::MatrixBase<MatrixType2>&>(Qaa)).noalias()
      += (dtau*dtau) * Pvv_next;
}

} // namespace idocp

#endif // IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_