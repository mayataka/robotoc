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
    const SplitRiccatiFactorization& riccati_next, const double dt, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  assert(dt > 0);
  if (has_floating_base_) {
    AtPqq_.template topRows<6>().noalias() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.Pqq.template topRows<6>();
    AtPqq_.bottomRows(dimv_-6) = riccati_next.Pqq.bottomRows(dimv_-6);
    AtPqv_.template topRows<6>().noalias() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.Pqv.template topRows<6>();
    AtPqv_.bottomRows(dimv_-6) = riccati_next.Pqv.bottomRows(dimv_-6);
    AtPvq_.template topRows<6>().noalias() 
        = kkt_matrix.Fqv().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.Pqq.template topRows<6>();
    AtPvq_.bottomRows(dimv_-6) = dt * riccati_next.Pqq.bottomRows(dimv_-6);
    AtPvv_.template topRows<6>().noalias() 
        = kkt_matrix.Fqv().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.Pqv.template topRows<6>();
    AtPvv_.bottomRows(dimv_-6) = dt * riccati_next.Pqv.bottomRows(dimv_-6);
  }
  else {
    AtPqq_ = riccati_next.Pqq;
    AtPqv_ = riccati_next.Pqv;
    AtPvq_ = dt * riccati_next.Pqq;
    AtPvv_ = dt * riccati_next.Pqv;
  }
  AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
  AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
  AtPvq_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pqv.transpose();
  AtPvv_.noalias() += kkt_matrix.Fvv().transpose() * riccati_next.Pvv;

  BtPq_.noalias() = kkt_matrix.Fvu.transpose() * riccati_next.Pqv.transpose();
  BtPv_.noalias() = kkt_matrix.Fvu.transpose() * riccati_next.Pvv;
  // Factorize F
  if (has_floating_base_) {
    kkt_matrix.Qqq().template leftCols<6>().noalias() 
        += AtPqq_.template leftCols<6>() 
            * kkt_matrix.Fqq().template topLeftCorner<6, 6>();
    kkt_matrix.Qqq().rightCols(dimv_-6).noalias() += AtPqq_.rightCols(dimv_-6);
    kkt_matrix.Qqv().template leftCols<6>().noalias() 
        += AtPqq_.template leftCols<6>() 
            * kkt_matrix.Fqv().template topLeftCorner<6, 6>();
    kkt_matrix.Qqv().rightCols(dimv_-6).noalias() 
        += dt * AtPqq_.rightCols(dimv_-6);
    kkt_matrix.Qvv().template leftCols<6>().noalias() 
        += AtPvq_.template leftCols<6>() 
            * kkt_matrix.Fqv().template topLeftCorner<6, 6>();
    kkt_matrix.Qvv().rightCols(dimv_-6).noalias() 
        += dt * AtPvq_.rightCols(dimv_-6);
  }
  else {
    kkt_matrix.Qqq().noalias() += AtPqq_;
    kkt_matrix.Qqv().noalias() += dt * AtPqq_;
    kkt_matrix.Qvv().noalias() += dt * AtPvq_;
  }
  kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
  kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
  kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  // Factorize H
  kkt_matrix.Qqu().noalias() += AtPqv_ * kkt_matrix.Fvu;
  kkt_matrix.Qvu().noalias() += AtPvv_ * kkt_matrix.Fvu;
  // Factorize G
  kkt_matrix.Quu.noalias() += BtPv_ * kkt_matrix.Fvu;
  // Factorize vector term
  kkt_residual.lu.noalias() += BtPq_ * kkt_residual.Fq();
  kkt_residual.lu.noalias() += BtPv_ * kkt_residual.Fv();
  kkt_residual.lu.noalias() -= kkt_matrix.Fvu.transpose() * riccati_next.sv;
}


inline void BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization(
    const SplitRiccatiFactorization& riccati_next, 
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const LQRPolicy& lqr_policy, const double dt, 
    SplitRiccatiFactorization& riccati) {
  assert(dt > 0);
  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pvv = kkt_matrix.Qvv();
  GK_.noalias() = kkt_matrix.Quu * lqr_policy.K; 
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
    riccati.sq.template head<6>().noalias() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.sq.template head<6>();
    riccati.sq.tail(dimv_-6) = riccati_next.sq.tail(dimv_-6);
    riccati.sv.template head<6>().noalias() 
        = kkt_matrix.Fqv().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.sq.template head<6>();
    riccati.sv.tail(dimv_-6) = dt * riccati_next.sq.tail(dimv_-6);
  }
  else {
    riccati.sq.noalias() = riccati_next.sq;
    riccati.sv.noalias() = dt * riccati_next.sq;
  }
  riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
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