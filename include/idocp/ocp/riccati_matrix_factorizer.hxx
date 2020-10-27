#ifndef IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_
#define IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_

#include "Eigen/LU"

#include <assert.h>

namespace idocp {

inline RiccatiMatrixFactorizer::RiccatiMatrixFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    AtPqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    BtPq_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BtPv_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())) {
}


inline RiccatiMatrixFactorizer::RiccatiMatrixFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dimu_(0),
    AtPqq_(),
    AtPqv_(),
    AtPvq_(),
    AtPvv_(),
    BtPq_(),
    BtPv_() {
}


inline RiccatiMatrixFactorizer::~RiccatiMatrixFactorizer() {
}


inline void RiccatiMatrixFactorizer::factorizeMatrices(
    const RiccatiFactorization& riccati_next, const double dtau, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual) {
  assert(dtau > 0);
  if (has_floating_base_) {
    AtPqq_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvq;
    AtPqv_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
    AtPvq_.noalias() = dtau * riccati_next.Pqq;
    AtPvq_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvq;
    AtPvv_.noalias() = dtau * riccati_next.Pqv;
    AtPvv_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  }
  else {
    AtPqq_.noalias() = riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvq;
    AtPqv_.noalias() = riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
    AtPvq_.noalias() = dtau * riccati_next.Pqq;
    AtPvq_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvq;
    AtPvv_.noalias() = dtau * riccati_next.Pqv;
    AtPvv_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  }
  BtPq_.noalias() = kkt_matrix.Fvu().transpose() * riccati_next.Pvq;
  BtPv_.noalias() = kkt_matrix.Fvu().transpose() * riccati_next.Pvv;
  // Factorize F
  if (has_floating_base_) {
    kkt_matrix.Qqq().noalias() += AtPqq_ * kkt_matrix.Fqq();
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
    kkt_matrix.Qqv().noalias() += dtau * AtPqq_;
    kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
    kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
    kkt_matrix.Qvv().noalias() += dtau * AtPvq_;
    kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
  }
  else {
    kkt_matrix.Qqq().noalias() += AtPqq_;
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
    kkt_matrix.Qqv().noalias() += dtau * AtPqq_;
    kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
    kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
    kkt_matrix.Qvv().noalias() += dtau * AtPvq_;
    kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
  }
  // Factorize H
  kkt_matrix.Qqu().noalias() += AtPqv_ * kkt_matrix.Fvu();
  kkt_matrix.Qvu().noalias() += AtPvv_ * kkt_matrix.Fvu();
  // Factorize G
  kkt_matrix.Quu().noalias() += BtPv_ * kkt_matrix.Fvu();
  // Factorize vector term
  kkt_residual.lu().noalias() += BtPq_ * kkt_residual.Fq();
  kkt_residual.lu().noalias() += BtPv_ * kkt_residual.Fv();
  kkt_residual.lu().noalias() -= kkt_matrix.Fvu().transpose() * riccati_next.sv;
}


inline void RiccatiMatrixFactorizer::factorizeRecursion(
    const RiccatiFactorization& riccati_next, const KKTMatrix& kkt_matrix, 
    const KKTResidual& kkt_residual, const RiccatiGain& gain, const double dtau, 
    RiccatiFactorization& riccati) const {
  assert(dtau > 0);
  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pvv = kkt_matrix.Qvv();
  riccati.Pqq.noalias() += kkt_matrix.Qqu() * gain.Kq();
  riccati.Pqv.noalias() += kkt_matrix.Qqu() * gain.Kv();
  riccati.Pvv.noalias() += kkt_matrix.Qvu() * gain.Kv();
  riccati.Pvq = riccati.Pqv.transpose();
  if (has_floating_base_) {
    riccati.sq.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
    riccati.sv.noalias() = dtau * riccati_next.sq;
    riccati.sv.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.sv;
  }
  else {
    riccati.sq.noalias() = riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
    riccati.sv.noalias() = dtau * riccati_next.sq;
    riccati.sv.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.sv;
  }
  riccati.sq.noalias() -= AtPqq_ * kkt_residual.Fq();
  riccati.sq.noalias() -= AtPqv_ * kkt_residual.Fv();
  riccati.sv.noalias() -= AtPvq_ * kkt_residual.Fq();
  riccati.sv.noalias() -= AtPvv_ * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();
  riccati.sv.noalias() -= kkt_residual.lv();
  riccati.sq.noalias() -= kkt_matrix.Qqu() * gain.k;
  riccati.sv.noalias() -= kkt_matrix.Qvu() * gain.k;
}

} // namespace idocp

#endif // IDOCP_RICCATI_MATRIX_FACTORIZER_HXX_