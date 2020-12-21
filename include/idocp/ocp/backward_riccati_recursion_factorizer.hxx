#ifndef IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 
#define IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/ocp/backward_riccati_recursion_factorizer.hpp"

#include <cassert>

namespace idocp {

inline BackwardRiccatiRecursionFactorizer::BackwardRiccatiRecursionFactorizer(
    const Robot& robot) 
  : has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    AtPqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    BtPq_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    BtPv_(Eigen::MatrixXd::Zero(robot.dimu(), robot.dimv())),
    GK_(Eigen::MatrixXd::Zero(robot.dimu(), 2*robot.dimv())) {
}


inline BackwardRiccatiRecursionFactorizer::BackwardRiccatiRecursionFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    dimu_(0),
    AtPqq_(),
    AtPqv_(),
    AtPvq_(),
    AtPvv_(),
    BtPq_(),
    BtPv_(),
    GK_() {
}


inline BackwardRiccatiRecursionFactorizer::
~BackwardRiccatiRecursionFactorizer() {
}


inline void BackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, const double dtau, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  assert(dtau >= 0);
  if (has_floating_base_) {
    AtPqq_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
    AtPqv_.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
  }
  else {
    AtPqq_ = riccati_next.Pqq;
    AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
    AtPqv_ = riccati_next.Pqv;
    AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
  }
  AtPvq_ = dtau * riccati_next.Pqq;
  AtPvq_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pqv.transpose();
  AtPvv_ = dtau * riccati_next.Pqv;
  AtPvv_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  BtPq_.noalias() = kkt_matrix.Fvu().transpose() * riccati_next.Pqv.transpose();
  BtPv_.noalias() = kkt_matrix.Fvu().transpose() * riccati_next.Pvv;
  // Factorize F
  if (has_floating_base_) {
    kkt_matrix.Qqq().noalias() += AtPqq_ * kkt_matrix.Fqq();
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
  }
  else {
    kkt_matrix.Qqq().noalias() += AtPqq_;
    kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
  }
  kkt_matrix.Qqv().noalias() += dtau * AtPqq_;
  kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qvv().noalias() += dtau * AtPvq_;
  kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
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


inline void BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization(
    const SplitRiccatiFactorization& riccati_next, 
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const LQRStateFeedbackPolicy& lqr_policy, const double dtau, 
    SplitRiccatiFactorization& riccati) {
  assert(dtau >= 0);
  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pvv = kkt_matrix.Qvv();
  GK_.noalias() = kkt_matrix.Quu() * lqr_policy.K; 
  riccati.Pqq.noalias() 
      -= lqr_policy.K.leftCols(dimv_).transpose() * GK_.leftCols(dimv_);
  riccati.Pqv.noalias() 
      -= lqr_policy.K.leftCols(dimv_).transpose() * GK_.rightCols(dimv_);
  riccati.Pvv.noalias() 
      -= lqr_policy.K.rightCols(dimv_).transpose() * GK_.rightCols(dimv_);
  riccati.Pvq = riccati.Pqv.transpose();
  // preserve the symmetry
  riccati.Pqq = 0.5 * (riccati.Pqq + riccati.Pqq.transpose()).eval();
  riccati.Pvv = 0.5 * (riccati.Pvv + riccati.Pvv.transpose()).eval();
  if (has_floating_base_) {
    riccati.sq.noalias() = kkt_matrix.Fqq().transpose() * riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
  }
  else {
    riccati.sq.noalias() = riccati_next.sq;
    riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
  }
  riccati.sv.noalias() = dtau * riccati_next.sq;
  riccati.sv.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.sv;
  riccati.sq.noalias() -= AtPqq_ * kkt_residual.Fq();
  riccati.sq.noalias() -= AtPqv_ * kkt_residual.Fv();
  riccati.sv.noalias() -= AtPvq_ * kkt_residual.Fq();
  riccati.sv.noalias() -= AtPvv_ * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();
  riccati.sv.noalias() -= kkt_residual.lv();
  riccati.sq.noalias() -= kkt_matrix.Qqu() * lqr_policy.k;
  riccati.sv.noalias() -= kkt_matrix.Qvu() * lqr_policy.k;
}

} // namespace idocp

#endif // IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 