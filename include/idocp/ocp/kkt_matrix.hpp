#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include <assert.h>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"


namespace idocp {

class KKTMatrix {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTMatrix(const Robot& robot) 
    : Quu(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      Fqq(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
      C_(Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), 
                               3*robot.dimv()+robot.max_dimf())),
      Q_(Eigen::MatrixXd::Zero(3*robot.dimv()+robot.max_dimf(), 
                               3*robot.dimv()+robot.max_dimf())),
      Sc_(Eigen::MatrixXd::Zero(robot.dim_passive()+robot.max_dimf(), 
                                robot.dim_passive()+robot.max_dimf())),
      Minv_(Eigen::MatrixXd::Zero(3*robot.dimv()+robot.dim_passive()+2*robot.max_dimf(), 
                                  3*robot.dimv()+robot.dim_passive()+2*robot.max_dimf())), 
      Sx_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())),
      FMinv_(Eigen::MatrixXd::Zero(2*robot.dimv(), 
                                   3*robot.dimv()+robot.dim_passive()+2*robot.max_dimf())),
      has_floating_base_(robot.has_floating_base()),
      dimv_(robot.dimv()), 
      dimx_(2*robot.dimv()), 
      dim_passive_(robot.dim_passive()), 
      max_dimf_(robot.max_dimf()), 
      dimf_(robot.dimf()), 
      max_dimc_(robot.dim_passive()+robot.max_dimf()), 
      dimc_(robot.dim_passive()+robot.dimf()),
      a_begin_(0),
      f_begin_(robot.dimv()),
      q_begin_(robot.dimv()+robot.dimf()),
      v_begin_(2*robot.dimv()+robot.dimf()),
      dimQ_(3*robot.dimv()+robot.dimf()) {
  }

  KKTMatrix() 
    : Quu(),
      Fqq(),
      C_(), 
      Q_(), 
      Sc_(), 
      Minv_(), 
      Sx_(), 
      FMinv_(),
      has_floating_base_(false),
      dimv_(0), 
      dimx_(0), 
      dim_passive_(0), 
      max_dimf_(0), 
      dimf_(0), 
      max_dimc_(0), 
      dimc_(0),
      a_begin_(0),
      f_begin_(0),
      q_begin_(0),
      v_begin_(0),
      dimQ_(0) {
  }

  ~KKTMatrix() {
  }

  KKTMatrix(const KKTMatrix&) = default;

  KKTMatrix& operator=(const KKTMatrix&) = default;
 
  KKTMatrix(KKTMatrix&&) noexcept = default;

  KKTMatrix& operator=(KKTMatrix&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    dimf_ = robot.dimf();
    dimc_ = robot.dim_passive() + robot.dimf();
    f_begin_ = robot.dimv();
    q_begin_ = robot.dimv() + robot.dimf();
    v_begin_ = 2*robot.dimv() + robot.dimf();
    dimQ_ = 3*robot.dimv() + robot.dimf();
  }

  inline Eigen::Block<Eigen::MatrixXd> Ca() {
    return C_.block(0, a_begin_, dimc_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Cf() {
    return C_.block(0, f_begin_, dimc_, dimf_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Cq() {
    return C_.block(0, q_begin_, dimc_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Cv() {
    return C_.block(0, v_begin_, dimc_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qaa() {
    return Q_.block(a_begin_, a_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qaf() {
    return Q_.block(a_begin_, f_begin_, dimv_, dimf_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qaq() {
    return Q_.block(a_begin_, q_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qav() {
    return Q_.block(a_begin_, v_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qfa() {
    return Q_.block(f_begin_, a_begin_, dimf_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qff() {
    return Q_.block(f_begin_, f_begin_, dimf_, dimf_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qfq() {
    return Q_.block(f_begin_, q_begin_, dimf_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qfv() {
    return Q_.block(f_begin_, v_begin_, dimf_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqa() {
    return Q_.block(q_begin_, a_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqf() {
    return Q_.block(q_begin_, f_begin_, dimv_, dimf_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqq() {
    return Q_.block(q_begin_, q_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qqv() {
    return Q_.block(q_begin_, v_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qva() {
    return Q_.block(v_begin_, a_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qvf() {
    return Q_.block(v_begin_, f_begin_, dimv_, dimf_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qvq() {
    return Q_.block(v_begin_, q_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qvv() {
    return Q_.block(v_begin_, v_begin_, dimv_, dimv_);
  }

  inline Eigen::Block<Eigen::MatrixXd> Qxx() {
    return Q_.block(q_begin_, q_begin_, dimx_, dimx_);
  }

  inline Eigen::Block<Eigen::MatrixXd> costHessian() {
    return Q_.topLeftCorner(dimQ_, dimQ_);
  }

  inline Eigen::Block<Eigen::MatrixXd> constraintsJacobian() {
    return C_.topLeftCorner(dimc_, dimQ_);
  }

  inline void symmetrize() {
    Q_.topLeftCorner(dimQ_, dimQ_).triangularView<Eigen::StrictlyLower>() 
        = Q_.topLeftCorner(dimQ_, dimQ_).transpose()
            .triangularView<Eigen::StrictlyLower>();
  }

  template <typename MatrixType>
  inline void invert(const double dtau,
                     const Eigen::MatrixBase<MatrixType>& kkt_matrix_inverse) {
    assert(kkt_matrix_inverse.rows() == (dimx_+dimc_+dimQ_));
    assert(kkt_matrix_inverse.cols() == (dimx_+dimc_+dimQ_));
    // Forms the Schur complement matrix
    invertConstrainedHessian(
        const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
            .block(dimx_, dimx_, dimQ_+dimc_, dimQ_+dimc_));
    if (has_floating_base_) {
      Sx_.topLeftCorner(dimv_, dimv_) 
          = Fqq * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, 
                                           dimx_+dimc_+q_begin_, dimv_, dimv_) 
                * Fqq.transpose()
            + dtau * Fqq 
                   * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, 
                                              dimx_+dimc_+v_begin_, dimv_, dimv_) 
            + dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                              dimx_+dimc_+q_begin_, dimv_, dimv_) 
                   * Fqq.transpose()
            + dtau * dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                                     dimx_+dimc_+v_begin_, dimv_, dimv_);
      Sx_.topRightCorner(dimv_, dimv_) 
          = - Fqq * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, 
                                             dimx_+dimc_+v_begin_, dimv_, dimv_)
            + dtau * Fqq 
                   * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, 
                                              dimx_+dimc_+a_begin_, dimv_, dimv_) 
            - dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                              dimx_+dimc_+v_begin_, dimv_, dimv_) 
            + dtau * dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                                     dimx_+dimc_+a_begin_, dimv_, dimv_);
    }
    else {
      Sx_.topLeftCorner(dimv_, dimv_) 
          = kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, dimx_+dimc_+q_begin_, 
                                     dimv_, dimv_)
              - dtau * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, 
                                                dimx_+dimc_+v_begin_, dimv_, dimv_) 
              - dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                                dimx_+dimc_+q_begin_, dimv_, dimv_)
              + dtau * dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_+dimc_+v_begin_, dimv_, dimv_);
      Sx_.topRightCorner(dimv_, dimv_) 
          = kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, dimx_+dimc_+v_begin_, 
                                     dimv_, dimv_)
              - dtau * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, 
                                                dimx_+dimc_+a_begin_, dimv_, dimv_)
              - dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                                dimx_+dimc_+v_begin_, dimv_, dimv_) 
              + dtau * dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                                       dimx_+dimc_+a_begin_, dimv_, dimv_);
    }
    Sx_.bottomLeftCorner(dimv_, dimv_) 
        = Sx_.topRightCorner(dimv_, dimv_).transpose();
    Sx_.bottomRightCorner(dimv_, dimv_) 
        = kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_+dimc_+v_begin_, 
                                   dimv_, dimv_)
            - dtau * kkt_matrix_inverse.block(dimx_+dimc_+a_begin_, 
                                              dimx_+dimc_+v_begin_, dimv_, dimv_) 
            - dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, 
                                              dimx_+dimc_+a_begin_, dimv_, dimv_) 
            + dtau * dtau * kkt_matrix_inverse.block(dimx_+dimc_+a_begin_, 
                                                     dimx_+dimc_+a_begin_, dimv_, dimv_);
    const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
        .topLeftCorner(dimx_, dimx_)
        = - Sx_.llt().solve(Eigen::MatrixXd::Identity(dimx_, dimx_));
    if (has_floating_base_) {
      FMinv_.topLeftCorner(dimv_, dimc_+dimQ_) 
          = dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_, 
                                            dimv_, dimc_+dimQ_)
              + Fqq * kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, dimx_, 
                                               dimv_, dimc_+dimQ_);
      FMinv_.bottomLeftCorner(dimv_, dimc_+dimQ_) 
          = dtau * kkt_matrix_inverse.block(dimx_+dimc_+a_begin_, dimx_, 
                                            dimv_, dimc_+dimQ_)
              - kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_, 
                                         dimv_, dimc_+dimQ_);
    }
    else {
      FMinv_.topLeftCorner(dimv_, dimc_+dimQ_) 
          = dtau * kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_, 
                                            dimv_, dimc_+dimQ_)
              - kkt_matrix_inverse.block(dimx_+dimc_+q_begin_, dimx_, 
                                          dimv_, dimc_+dimQ_);
      FMinv_.bottomLeftCorner(dimv_, dimc_+dimQ_) 
          = dtau * kkt_matrix_inverse.block(dimx_+dimc_+a_begin_, dimx_, 
                                            dimv_, dimc_+dimQ_)
              - kkt_matrix_inverse.block(dimx_+dimc_+v_begin_, dimx_, 
                                         dimv_, dimc_+dimQ_);
    }
    const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
        .block(0, dimx_, dimx_, dimc_+dimQ_)
        = - kkt_matrix_inverse.topLeftCorner(dimx_, dimx_)
            * FMinv_.topLeftCorner(dimx_, dimc_+dimQ_);
    const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
        .block(dimx_, 0, dimc_+dimQ_, dimx_)
        = kkt_matrix_inverse.block(0, dimx_, dimx_, dimc_+dimQ_).transpose();
    const_cast<Eigen::MatrixBase<MatrixType>&>(kkt_matrix_inverse)
        .block(dimx_, dimx_, dimc_+dimQ_, dimc_+dimQ_)
        -= kkt_matrix_inverse.block(0, dimx_, dimx_, dimc_+dimQ_).transpose()
                * Sx_ * kkt_matrix_inverse.block(0, dimx_, dimx_, dimc_+dimQ_);
  }

  inline void setZeroMinimum() {
    Quu.setZero();
    Fqq.topLeftCorner<6, 6>().setZero();
    C_.topLeftCorner(dimc_, dimQ_).setZero();
    Q_.topLeftCorner(dimQ_, dimQ_).triangularView<Eigen::Upper>().setZero();
  }

  inline void setZero() {
    Quu.setZero();
    Fqq.setZero();
    C_.setZero();
    Q_.setZero();
  }

  Eigen::MatrixXd Quu, Fqq;

private:
  Eigen::MatrixXd C_, Q_, Sc_, Minv_, Sx_, FMinv_;
  bool has_floating_base_;
  int dimv_, dimx_, dim_passive_, max_dimf_, dimf_, max_dimc_, dimc_, 
      a_begin_, f_begin_, q_begin_, v_begin_, dimQ_;

  template <typename MatrixType>
  inline void invertConstrainedHessian( 
      const Eigen::MatrixBase<MatrixType>& hessian_inverse) {
    assert(hessian_inverse.rows() == dimQ_+dimc_);
    assert(hessian_inverse.cols() == dimQ_+dimc_);
    if (dimc_ > 0) {
      const_cast<Eigen::MatrixBase<MatrixType>&>(hessian_inverse)
          .block(dimc_, dimc_, dimQ_, dimQ_)
          = Q_.topLeftCorner(dimQ_, dimQ_)
              .llt().solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
      Sc_.topLeftCorner(dimc_, dimc_) 
          = C_.topLeftCorner(dimc_, dimQ_) 
              * hessian_inverse.block(dimc_, dimc_, dimQ_, dimQ_)
              * C_.topLeftCorner(dimc_, dimQ_).transpose();
      const_cast<Eigen::MatrixBase<MatrixType>&>(hessian_inverse)
          .topLeftCorner(dimc_, dimc_)
          = Sc_.topLeftCorner(dimc_, dimc_)
               .llt().solve(Eigen::MatrixXd::Identity(dimc_, dimc_));
      const_cast<Eigen::MatrixBase<MatrixType>&>(hessian_inverse)
          .block(0, dimc_, dimc_, dimQ_)
          = hessian_inverse.topLeftCorner(dimc_, dimc_) 
              * C_.topLeftCorner(dimc_, dimQ_) 
              * hessian_inverse.block(dimc_, dimc_, dimQ_, dimQ_);
      const_cast<Eigen::MatrixBase<MatrixType>&>(hessian_inverse)
          .block(dimc_, 0, dimQ_, dimc_)
          = hessian_inverse.block(0, dimc_, dimc_, dimQ_).transpose();
      const_cast<Eigen::MatrixBase<MatrixType>&>(hessian_inverse)
          .block(dimc_, dimc_, dimQ_, dimQ_).noalias()
          -= hessian_inverse.block(0, dimc_, dimc_, dimQ_).transpose()
                * Sc_.topLeftCorner(dimc_, dimc_)
                * hessian_inverse.block(0, dimc_, dimc_, dimQ_);
    }
    else {
      const_cast<Eigen::MatrixBase<MatrixType>&>(hessian_inverse)
          = Q_.topLeftCorner(dimQ_, dimQ_)
              .llt().solve(Eigen::MatrixXd::Identity(dimQ_, dimQ_));
    }
  }

};

} // namespace idocp 


#endif // IDOCP_KKT_MATRIX_HPP_