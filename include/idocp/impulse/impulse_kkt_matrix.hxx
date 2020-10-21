#ifndef IDOCP_IMPULSE_KKT_MATRIX_HXX_
#define IDOCP_IMPULSE_KKT_MATRIX_HXX_

#include "idocp/impulse/impulse_kkt_matrix.hpp"

#include <assert.h>

#include "Eigen/LU"

namespace idocp {

inline ImpulseKKTMatrix::ImpulseKKTMatrix(const Robot& robot)
  : Qdvdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fvq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fvv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq_prev(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    schur_complement_(2*robot.max_dimf(), 2*robot.dimv()+robot.max_dimf()),
    C_(Eigen::MatrixXd::Zero(robot.max_dimf(), 
                             2*robot.dimv()+robot.max_dimf())),
    Q_(Eigen::MatrixXd::Zero(2*robot.dimv()+robot.max_dimf(), 
                             2*robot.dimv()+robot.max_dimf())),
    Fvf_full_(Eigen::MatrixXd::Zero(robot.dimv(), robot.max_dimf())),
    Sx_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    FMinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 
                                 2*robot.dimv()+2*robot.max_dimf())),
    has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimc_(0),
    dimKKT_(4*robot.dimv()),
    dimQ_(2*robot.dimv()),
    f_begin_(0),
    q_begin_(0),
    v_begin_(robot.dimv()),
    max_dimKKT_(4*robot.dimv()+2*robot.max_dimf()) {
}
 

inline ImpulseKKTMatrix::ImpulseKKTMatrix() 
  : Qdvdv(),
    Fqq(),
    Fvq(),
    Fvv(),
    Fqq_prev(),
    schur_complement_(),
    C_(), 
    Q_(), 
    Fvf_full_(),
    Sx_(), 
    FMinv_(),
    has_floating_base_(false),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimc_(0),
    dimKKT_(0),
    dimQ_(0),
    f_begin_(0),
    q_begin_(0),
    v_begin_(0),
    max_dimKKT_(0) {
}


inline ImpulseKKTMatrix::~ImpulseKKTMatrix() {
}


inline void ImpulseKKTMatrix::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  dimc_ = contact_status.dimf();
  dimKKT_ = 4*dimv_ + dimf_ + dimc_;
  dimQ_ = 2*dimv_ + dimf_;
  q_begin_ = contact_status.dimf();
  v_begin_ = contact_status.dimf() + dimv_;
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Cf() {
  return C_.block(0, f_begin_, dimc_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Cq() {
  return C_.block(0, q_begin_, dimc_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Cv() {
  return C_.block(0, v_begin_, dimc_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Cqv() {
  return C_.block(0, q_begin_, dimc_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qff() {
  return Q_.block(f_begin_, f_begin_, dimf_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qfq() {
  return Q_.block(f_begin_, q_begin_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qfv() {
  return Q_.block(f_begin_, v_begin_, dimf_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqf() {
  return Q_.block(q_begin_, f_begin_, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqq() {
  return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqv() {
  return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qvf() {
  return Q_.block(v_begin_, f_begin_, dimv_, dimf_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qvq() {
  return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qvv() {
  return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qxx() {
  return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::costHessian() {
  return Q_.topLeftCorner(dimQ_, dimQ_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::constraintsJacobian() {
  return C_.topLeftCorner(dimc_, dimQ_);
}


inline void ImpulseKKTMatrix::symmetrize() {
  Q_.topLeftCorner(dimQ_, dimQ_).triangularView<Eigen::StrictlyLower>() 
    = Q_.topLeftCorner(dimQ_, dimQ_).transpose()
                                    .triangularView<Eigen::StrictlyLower>();
}


template <typename MatrixType>
inline void ImpulseKKTMatrix::invert(
    const Eigen::MatrixBase<MatrixType>& kkt_matrix_inverse) {
  // Forms the Schur complement matrix
  const int dimcQ = dimc_ + dimQ_;
  assert(kkt_matrix_inverse.rows() == (dimx_+dimcQ));
  assert(kkt_matrix_inverse.cols() == (dimx_+dimcQ));
  invertConstrainedHessian(
      const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
          .bottomRightCorner(dimcQ, dimcQ));
  FMinv_.setZero();
  if (has_floating_base_) {
    FMinv_.topLeftCorner(dimv_, dimcQ).template topRows<kDimFloatingBase>()
        = Fqq.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>() 
            * kkt_matrix_inverse.block(dimx_+dimc_+dimf_, dimx_, kDimFloatingBase, dimcQ);
    FMinv_.topLeftCorner(dimv_, dimcQ).bottomRows(dimv_-kDimFloatingBase)
        = - kkt_matrix_inverse.block(dimx_+dimc_+dimf_+kDimFloatingBase, dimx_,
                                     dimv_-kDimFloatingBase, dimcQ);
  }
  else {
    FMinv_.topLeftCorner(dimv_, dimcQ) 
        = - kkt_matrix_inverse.block(dimx_+dimc_+dimf_, dimx_, dimv_, dimcQ);
  }
  FMinv_.bottomLeftCorner(dimv_, dimcQ).noalias()
      = Fvf() * kkt_matrix_inverse.block(dimx_+dimc_, dimx_, dimf_, dimcQ);
  FMinv_.bottomLeftCorner(dimv_, dimcQ).noalias()
      += Fvq * kkt_matrix_inverse.block(dimx_+dimc_+dimf_, dimx_, dimv_, dimcQ);
  FMinv_.bottomLeftCorner(dimv_, dimcQ).noalias()
      -= kkt_matrix_inverse.block(dimx_+dimc_+dimf_+dimv_, dimx_, dimv_, dimcQ);
  std::cout << "FMinv_" << std::endl;
  std::cout << FMinv_ << std::endl;
  if (has_floating_base_) {
    Sx_.topLeftCorner(dimv_, dimv_).template leftCols<kDimFloatingBase>().noalias()
        += FMinv_.block(0, dimc_+dimf_, dimv_, dimv_).template leftCols<kDimFloatingBase>() 
            * Fqq.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>().transpose();
    Sx_.topLeftCorner(dimv_, dimv_).rightCols(dimv_-kDimFloatingBase).noalias()
        -= FMinv_.block(0, dimc_+dimf_, dimv_, dimv_).rightCols(dimv_-kDimFloatingBase);
    Sx_.bottomLeftCorner(dimv_, dimv_).template leftCols<kDimFloatingBase>().noalias()
        += FMinv_.block(dimv_, dimc_+dimf_, dimv_, dimv_).template leftCols<kDimFloatingBase>() 
            * Fqq.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>().transpose();
    Sx_.bottomLeftCorner(dimv_, dimv_).rightCols(dimv_-kDimFloatingBase).noalias()
        -= FMinv_.block(dimv_, dimc_+dimf_, dimv_, dimv_).rightCols(dimv_-kDimFloatingBase);
  }
  else {
    Sx_.leftCols(dimv_) = - FMinv_.block(0, dimc_+dimf_, dimx_, dimv_);
  }
  Sx_.rightCols(dimv_).noalias() = FMinv_.block(0, dimc_, dimx_, dimf_) * Fvf().transpose();
  Sx_.rightCols(dimv_).noalias() += FMinv_.block(0, dimc_+dimf_, dimx_, dimv_) * Fvq.transpose();
  Sx_.rightCols(dimv_).noalias() -= FMinv_.block(0, dimc_+dimf_+dimv_, dimx_, dimv_);
  std::cout << "Sx_" << std::endl;
  std::cout << Sx_ << std::endl;
  const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
      .topLeftCorner(dimx_, dimx_).noalias()
      = - Sx_.llt().solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
      .topRightCorner(dimx_, dimcQ).noalias()
      = - kkt_matrix_inverse.topLeftCorner(dimx_, dimx_)
          * FMinv_.topLeftCorner(dimx_, dimcQ);
  const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
      .bottomLeftCorner(dimcQ, dimx_)
      = kkt_matrix_inverse.topRightCorner(dimx_, dimcQ).transpose();
  const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
      .bottomRightCorner(dimcQ, dimcQ).noalias()
      -= kkt_matrix_inverse.topRightCorner(dimx_, dimcQ).transpose()
              * Sx_ * kkt_matrix_inverse.topRightCorner(dimx_, dimcQ);
}


inline void ImpulseKKTMatrix::setZero() {
  Qdvdv.setZero();
  Fqq.setZero();
  Fvq.setZero();
  Fvv.setZero();
  Fqq_prev.setZero();
  C_.setZero();
  Q_.setZero();
}


inline int ImpulseKKTMatrix::dimKKT() const {
  return dimKKT_;
}


inline int ImpulseKKTMatrix::max_dimKKT() const {
  return max_dimKKT_;
}


inline int ImpulseKKTMatrix::dimc() const {
  return dimc_;
}


inline int ImpulseKKTMatrix::dimf() const {
  return dimf_;
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Fvf() {
  return Fvf_full_.block(0, 0, dimv_, dimf_);
}


template <typename MatrixType>
inline void ImpulseKKTMatrix::invertConstrainedHessian( 
    const Eigen::MatrixBase<MatrixType>& Hinv) {
  assert(dimc_ > 0);
  schur_complement_.invertWithZeroTopLeftCorner(
      dimc_, dimQ_, C_.topLeftCorner(dimc_, dimQ_),
      Q_.topLeftCorner(dimQ_, dimQ_), 
      const_cast<Eigen::MatrixBase<MatrixType>&>(Hinv));
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_KKT_MATRIX_HXX_ 