#ifndef IDOCP_RICCATI_FACTORIZER_HXX_
#define IDOCP_RICCATI_FACTORIZER_HXX_

#include "idocp/ocp/riccati_factorizer.hpp"

#include <assert.h>

namespace idocp {

inline RiccatiFactorizer::RiccatiFactorizer(const Robot& robot) 
  : has_floating_base_(robot.has_floating_base()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    llt_(robot.dimu()),
    AtPqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    BtPq_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BtPv_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    Ginv_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimu())) {
}


inline RiccatiFactorizer::RiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dimu_(0),
    llt_(),
    AtPqq_(),
    AtPqv_(),
    AtPvq_(),
    AtPvv_(),
    BtPq_(),
    BtPv_(),
    Ginv_() {
}


inline RiccatiFactorizer::~RiccatiFactorizer() {
}


inline void RiccatiFactorizer::factorizeBackwardRicursion(
    const RiccatiSolution& riccati_next, const double dtau, 
    KKTMatrix& kkt_matrix, KKTResidual& kkt_residual, RiccatiGain& gain, 
    RiccatiSolution& riccati) {
  assert(dtau > 0);
  factorizeMatrices(riccati_next, dtau, kkt_matrix, kkt_residual);
  computeFeedbackGainAndFeedforward(kkt_matrix, kkt_residual, gain);
  factorizeRecursion(riccati_next, dtau, kkt_matrix, kkt_residual, gain, riccati);
}


inline void RiccatiFactorizer::factorizeForwardRicursion(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual,
    const SplitDirection& d, const double dtau, SplitDirection& d_next) const {
  assert(dtau > 0);
  if (has_floating_base_) {
    d_next.dq().noalias() = kkt_matrix.Fqq() * d.dq() + dtau * d.dv() 
                              + kkt_residual.Fq();
  }
  else {
    d_next.dq().noalias() = d.dq() + dtau * d.dv() + kkt_residual.Fq();
  }
  d_next.dv().noalias() = kkt_matrix.Fvq() * d.dq() + kkt_matrix.Fvv() * d.dv() 
                            + kkt_matrix.Fvu() * d.du() + kkt_residual.Fv();
}


inline void RiccatiFactorizer::computeCostateDirection(
    const RiccatiSolution& riccati, SplitDirection& d) {
  d.dlmd().noalias() = riccati.Pqq * d.dq() + riccati.Pqv * d.dv() - riccati.sq;
  d.dgmm().noalias() = riccati.Pvq * d.dq() + riccati.Pvv * d.dv() - riccati.sv;
}


inline void RiccatiFactorizer::computeControlInputDirection(
    const RiccatiGain& gain, SplitDirection& d) {
  d.du() = gain.k;
  d.du().noalias() += gain.K * d.dx();
}


inline void RiccatiFactorizer::factorizeMatrices(
    const RiccatiSolution& riccati_next, const double dtau, 
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


inline void RiccatiFactorizer::computeFeedbackGainAndFeedforward(
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
    RiccatiGain& gain) {
  llt_.compute(kkt_matrix.Quu());
  Ginv_ = llt_.solve(Eigen::MatrixXd::Identity(dimu_, dimu_));
  gain.K.noalias() = - Ginv_ * kkt_matrix.Qxu().transpose();
  gain.k.noalias() = - Ginv_ * kkt_residual.lu();
}


inline void RiccatiFactorizer::factorizeRecursion(
    const RiccatiSolution& riccati_next, const double dtau, 
    const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
    const RiccatiGain& gain, RiccatiSolution& riccati) const {
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

#endif // IDOCP_RICCATI_FACTORIZER_HXX_