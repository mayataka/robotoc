#ifndef IDOCP_IMPULSE_RICCATI_MATRIX_FACTORIZER_HXX_
#define IDOCP_IMPULSE_RICCATI_MATRIX_FACTORIZER_HXX_

#include "idocp/impulse/impulse_riccati_matrix_factorizer.hpp"

namespace idocp {

inline ImpulseRiccatiMatrixFactorizer::ImpulseRiccatiMatrixFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    ATPqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    ATPqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    ATPvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    ATPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
}


inline ImpulseRiccatiMatrixFactorizer::ImpulseRiccatiMatrixFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    ATPqq_(),
    ATPqv_(),
    ATPvq_(),
    ATPvv_() {
}


inline ImpulseRiccatiMatrixFactorizer::~ImpulseRiccatiMatrixFactorizer() {
}


inline void ImpulseRiccatiMatrixFactorizer::factorize(
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual,
    const RiccatiFactorization& riccati_next, RiccatiFactorization& riccati) {
  ATPqq_.noalias() = kkt_matrix.Fqq.transpose() * riccati_next.Pqq;
  ATPqq_.noalias() += kkt_matrix.Fvq.transpose() * riccati_next.Pvq;
  ATPqv_.noalias() = kkt_matrix.Fqq.transpose() * riccati_next.Pqv;
  ATPqv_.noalias() += kkt_matrix.Fvq.transpose() * riccati_next.Pvv;
  ATPvq_.noalias() = kkt_matrix.Fvv.transpose() * riccati_next.Pvq;
  ATPvv_.noalias() = kkt_matrix.Fvv.transpose() * riccati_next.Pvv;

  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqq.noalias() += ATPqq_ * kkt_matrix.Fqq;
  riccati.Pqq.noalias() += ATPqv_ * kkt_matrix.Fvq;

  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pqv.noalias() += ATPqv_ * kkt_matrix.Fvv;

  riccati.Pvq = riccati.Pqv.transpose();

  riccati.Pvv = kkt_matrix.Qvv();
  riccati.Pvv.noalias() += ATPvv_ * kkt_matrix.Fvv;

  riccati.sq.noalias() = kkt_matrix.Fqq.transpose() * riccati_next.sq;
  riccati.sq.noalias() += kkt_matrix.Fvq.transpose() * riccati_next.sv;
  riccati.sq.noalias() -= ATPqq_ * kkt_residual.Fq();
  riccati.sq.noalias() -= ATPqv_ * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();

  riccati.sv.noalias() = kkt_matrix.Fvv.transpose() * riccati_next.sv;
  riccati.sv.noalias() -= ATPvq_ * kkt_residual.Fq();
  riccati.sv.noalias() -= ATPvv_ * kkt_residual.Fv();
  riccati.sv.noalias() -= kkt_residual.lv();
}

} // namespace idocp

#endif // IDOCP_IMPULSE_RICCATI_MATRIX_FACTORIZER_HXX_