#ifndef IDOCP_IMPULSE_KKT_MATRIX_HXX_
#define IDOCP_IMPULSE_KKT_MATRIX_HXX_

#include "idocp/impulse/impulse_kkt_matrix.hpp"

#include <assert.h>

#include "Eigen/LU"

namespace idocp {

inline ImpulseKKTMatrix::ImpulseKKTMatrix(
    const Robot& robot, const bool use_contact_position_constraint)
  : Qdvdv(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    Fqq_prev(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    C_(Eigen::MatrixXd::Zero(robot.max_dimf(), 2*robot.dimv())),
    Q_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    Sc_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.max_dimf())),
    Sx_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
    FMinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 
                                 2*robot.dimv()+robot.max_dimf())),
    C_H_inv_(Eigen::MatrixXd::Zero(robot.max_dimf(), 
                                   2*robot.dimv()+robot.max_dimf())),
    Fqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(4*robot.dimv()+3*robot.max_dimf()),
    use_contact_position_constraint_(use_contact_position_constraint) {
}
 

inline ImpulseKKTMatrix::ImpulseKKTMatrix() 
  : Qdvdv(),
    Fqq(),
    Fqq_prev(),
    C_(), 
    Q_(), 
    Sc_(), 
    Sx_(), 
    FMinv_(),
    C_H_inv_(),
    has_floating_base_(false),
    dimv_(0), 
    dimx_(0), 
    dimf_(0), 
    dimc_(0),
    max_dimKKT_(0),
    use_contact_position_constraint_(false) {
}


inline ImpulseKKTMatrix::~ImpulseKKTMatrix() {
}


inline void ImpulseKKTMatrix::setContactStatus(
    const ContactStatus& contact_status) {
  dimf_ = contact_status.dimf();
  dimc_ = 2*contact_status.dimf();
  dimKKT_ = 4*dimv_ + dimf_ + dimc_;
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


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqq() {
  return Q_.block(0, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qqv() {
  return Q_.block(0, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qvq() {
  return Q_.block(dimv_, 0, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qvv() {
  return Q_.block(dimv_, dimv_, dimv_, dimv_);
}


inline Eigen::Block<Eigen::MatrixXd> ImpulseKKTMatrix::Qxx() {
  return Q_;
}


inline void ImpulseKKTMatrix::symmetrize() {
  Q_.triangularView<Eigen::StrictlyLower>() 
    = Q_.transpose().triangularView<Eigen::StrictlyLower>();
}


template <typename MatrixType>
inline void ImpulseKKTMatrix::invert(
    const Eigen::MatrixBase<MatrixType>& kkt_matrix_inverse) {
  assert(kkt_matrix_inverse.rows() == (2*dimx_+dimc_));
  assert(kkt_matrix_inverse.cols() == (2*dimx_+dimc_));
  // Forms the Schur complement matrix
  const int dimcQ = dimc_ + dimx_;
  invertConstrainedHessian(
      const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
          .bottomRightCorner(dimcQ, dimcQ));
  if (has_floating_base_) {
    FMinv_.topLeftCorner(dimv_, dimcQ) 
        = dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_, dimv_, dimcQ);
    FMinv_.topLeftCorner(dimv_, dimcQ).template topRows<kDimFloatingBase>().noalias()
        += Fqq.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>() 
            * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, dimx_, 
                                       kDimFloatingBase, dimcQ);
    FMinv_.topLeftCorner(dimv_, dimcQ).bottomRows(dimv_-kDimFloatingBase).noalias()
        -= kkt_matrix_inverse.block(dimx_+dimc_+q_begin_+kDimFloatingBase, 
                                    dimx_, dimv_-kDimFloatingBase, dimcQ);
  }
  else {
    FMinv_.topLeftCorner(dimv_, dimcQ).noalias() 
        = dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_, 
                                          dimv_, dimcQ)
          - kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, dimx_, dimv_, dimcQ);
  }
  FMinv_.bottomLeftCorner(dimv_, dimcQ).noalias() 
      = dtau * kkt_matrix_inverse.block(dimx_+dimc_+f_begin_, dimx_, 
                                        dimv_, dimcQ)
        - kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_, dimv_, dimcQ);
  if (has_floating_base_) {
    Sx_.topLeftCorner(dimv_, dimv_) 
        = dtau * FMinv_.block(0, dimc_+v_begin_, dimv_, dimv_);
    Sx_.topLeftCorner(dimv_, dimv_).template leftCols<kDimFloatingBase>().noalias()
        += FMinv_.block(0, dimc_+q_begin_, 
                        dimv_, dimv_).template leftCols<kDimFloatingBase>() 
            * Fqq.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>().transpose();
    Sx_.topLeftCorner(dimv_, dimv_).rightCols(dimv_-kDimFloatingBase).noalias()
        -= FMinv_.block(0, dimc_+q_begin_, 
                        dimv_, dimv_).rightCols(dimv_-kDimFloatingBase);
    Sx_.bottomLeftCorner(dimv_, dimv_) 
        = dtau * FMinv_.block(dimv_, dimc_+v_begin_, dimv_, dimv_);
    Sx_.bottomLeftCorner(dimv_, dimv_).template leftCols<kDimFloatingBase>().noalias() 
        += FMinv_.block(dimv_, dimc_+q_begin_, 
                        dimv_, dimv_).template leftCols<kDimFloatingBase>() 
            * Fqq.template topLeftCorner<kDimFloatingBase, kDimFloatingBase>().transpose();
    Sx_.bottomLeftCorner(dimv_, dimv_).rightCols(dimv_-kDimFloatingBase).noalias() 
        -= FMinv_.block(dimv_, dimc_+q_begin_, 
                        dimv_, dimv_).rightCols(dimv_-kDimFloatingBase);
  }
  else {
    Sx_.topLeftCorner(dimv_, dimv_).noalias() 
        = dtau * FMinv_.block(0, dimc_+v_begin_, dimv_, dimv_) 
            - FMinv_.block(0, dimc_+q_begin_, dimv_, dimv_);
    Sx_.bottomLeftCorner(dimv_, dimv_).noalias() 
        = dtau * FMinv_.block(dimv_, dimc_+v_begin_, dimv_, dimv_) 
            - FMinv_.block(dimv_, dimc_+q_begin_, dimv_, dimv_);
  }
  Sx_.topRightCorner(dimv_, dimv_).noalias() 
      = dtau * FMinv_.block(0, dimc_+f_begin_, dimv_, dimv_)
        - FMinv_.block(0, dimc_+v_begin_, dimv_, dimv_);
  Sx_.bottomRightCorner(dimv_, dimv_).noalias() 
      = dtau * FMinv_.block(dimv_, dimc_+f_begin_, dimv_, dimv_)
        - FMinv_.block(dimv_, dimc_+v_begin_, dimv_, dimv_);
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
  C_.setZero();
  Q_.setZero();
  Qff_full_.setZero();
  C_contact_velocity_full_.setZero();
}


inline int ImpulseKKTMatrix::dimKKT() const {
  return 4*dimv_+dimf_+dimc_;
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


template <typename MatrixType>
inline void ImpulseKKTMatrix::invertConstrainedHessian( 
    const Eigen::MatrixBase<MatrixType>& H_inv) {
  assert(H_inv.rows() == dimx_+dimc_);
  assert(H_inv.cols() == dimx_+dimc_);
  const_cast<Eigen::MatrixBase<MatrixType>&>(H_inv)
      .bottomRightCorner(dimx_, dimx_).noalias()
      = Q_.llt().solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
  C_H_inv_.topLeftCorner(dimc_, dimx_).noalias()
      = C_.topLeftCorner(dimc_, dimx_) * H_inv.bottomRightCorner(dimx_, dimx_);
  Sc_.topLeftCorner(dimc_, dimc_).noalias() 
      = C_H_inv_.topLeftCorner(dimc_, dimx_)
          * C_.topLeftCorner(dimc_, dimx_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType>&>(H_inv)
      .topLeftCorner(dimc_, dimc_).noalias()
      = - Sc_.topLeftCorner(dimc_, dimc_)
              .llt().solve(Eigen::MatrixXd::Identity(dimc_, dimc_));
  const_cast<Eigen::MatrixBase<MatrixType>&>(H_inv)
      .topRightCorner(dimc_, dimx_).noalias()
      = - H_inv.topLeftCorner(dimc_, dimc_) 
          * C_H_inv_.topLeftCorner(dimc_, dimx_);
  const_cast<Eigen::MatrixBase<MatrixType>&>(H_inv)
      .bottomLeftCorner(dimx_, dimc_)
      = H_inv.topRightCorner(dimc_, dimx_).transpose();
  const_cast<Eigen::MatrixBase<MatrixType>&>(H_inv)
      .bottomRightCorner(dimx_, dimx_).noalias()
      -= H_inv.topRightCorner(dimc_, dimx_).transpose()
            * Sc_.topLeftCorner(dimc_, dimc_)
            * H_inv.topRightCorner(dimc_, dimx_);
}

} // namespace idocp 

#endif // IDOCP_IMPULSE_KKT_MATRIX_HXX_ 