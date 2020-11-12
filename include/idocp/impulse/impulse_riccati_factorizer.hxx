#ifndef IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_
#define IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_

#include "idocp/impulse/impulse_riccati_factorizer.hpp"

namespace idocp {

inline ImpulseRiccatiFactorizer::ImpulseRiccatiFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    AtPqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    NApBKt_(Eigen::MatrixXd::Zero(2*robot.dimv(), 2*robot.dimv())) {
}


inline ImpulseRiccatiFactorizer::ImpulseRiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    AtPqq_(),
    AtPqv_(),
    AtPvq_(),
    AtPvv_(),
    NApBKt_() {
}


inline ImpulseRiccatiFactorizer::~ImpulseRiccatiFactorizer() {
}


inline void ImpulseRiccatiFactorizer::backwardRiccatiRecursion(
    const RiccatiFactorization& riccati_next, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual, RiccatiFactorization& riccati) {
  factorizeKKTMatrix(riccati_next, kkt_matrix, kkt_residual);
  factorizeRiccatiFactorization(riccati_next, kkt_matrix, kkt_residual, 
                                riccati);
}


inline void ImpulseRiccatiFactorizer::forwardStateConstraintFactorizationSerial(
    const RiccatiFactorization& riccati, const ImpulseKKTMatrix& kkt_matrix, 
    const ImpulseKKTResidual& kkt_residual, 
    RiccatiFactorization& riccati_next) {
  riccati_next.Pi.noalias() = kkt_matrix.Fxx() * riccati.Pi;
  riccati_next.pi = kkt_residual.Fx();
  riccati_next.pi.noalias() += kkt_matrix.Fxx() * riccati.pi;
  NApBKt_.noalias() = riccati.N * kkt_matrix.Fxx().transpose();
  riccati_next.N.noalias() = kkt_matrix.Fxx() * NApBKt_;
}


template <typename SplitDirectionType>
inline void ImpulseRiccatiFactorizer::forwardRiccatiRecursionSerial(
    const RiccatiFactorization& riccati, const ImpulseKKTMatrix& kkt_matrix, 
    const ImpulseKKTResidual& kkt_residual, const ImpulseSplitDirection& d, 
    SplitDirectionType& d_next) const {
  assert(dtau > 0);
  if (has_floating_base_) {
    d_next.dq().noalias() = kkt_matrix.Fqq() * d.dq();
    d_next.dq().noalias() += kkt_residual.Fq();
  }
  else {
    d_next.dq().noalias() = d.dq() + kkt_residual.Fq();
  }
  d_next.dv().noalias() = kkt_matrix.Fvq() * d.dq();
  d_next.dv().noalias() += kkt_matrix.Fvv() * d.dv();
  d_next.dv().noalias() += kkt_residual.Fv();
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


inline void ImpulseRiccatiFactorizer::factorizeKKTMatrix(
    const RiccatiFactorization& riccati_next, ImpulseKKTMatrix& kkt_matrix, 
    ImpulseKKTResidual& kkt_residual) {
  if (has_floating_base_) {
    AtPqq_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
    AtPqv_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
    AtPvq_.noalias() = kkt_matrix.Fvv().transpose() * riccati_next.Pqv.transpose();
    AtPvv_.noalias() = kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  }
  else {
    AtPqq_.noalias() = riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
    AtPqv_.noalias() = riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
    AtPvq_.noalias() = kkt_matrix.Fvv().transpose() * riccati_next.Pqv.transpose();
    AtPvv_.noalias() = kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  }
  // Factorize F
  if (has_floating_base_) {
    kkt_matrix.Qqq().noalias() += AtPqq_ * kkt_matrix.Fqq();
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
    kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
    kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
    kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
  }
  else {
    kkt_matrix.Qqq().noalias() += AtPqq_;
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
    kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
    kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
    kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
  }
}


inline void ImpulseRiccatiFactorizer::factorizeRiccatiFactorization(
    const RiccatiFactorization& riccati_next, const ImpulseKKTMatrix& kkt_matrix, 
    const ImpulseKKTResidual& kkt_residual, RiccatiFactorization& riccati) {
  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pvv = kkt_matrix.Qvv();
  riccati.Pvq = riccati.Pqv.transpose();
  // preserve the symmetry
  riccati.Pqq = 0.5 * (riccati.Pqq + riccati.Pqq.transpose()).eval();
  riccati.Pvv = 0.5 * (riccati.Pvv + riccati.Pvv.transpose()).eval();
  if (has_floating_base_) {
    riccati.sq.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
    riccati.sv.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.sv;
  }
  else {
    riccati.sq.noalias() = riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
    riccati.sv.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.sv;
  }
  riccati.sq.noalias() -= AtPqq_ * kkt_residual.Fq();
  riccati.sq.noalias() -= AtPqv_ * kkt_residual.Fv();
  riccati.sv.noalias() -= AtPvq_ * kkt_residual.Fq();
  riccati.sv.noalias() -= AtPvv_ * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();
  riccati.sv.noalias() -= kkt_residual.lv();
}

} // namespace idocp

#endif // IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_