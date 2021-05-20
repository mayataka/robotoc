#ifndef IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_ 
#define IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_

#include "idocp/riccati/impulse_backward_riccati_recursion_factorizer.hpp"

namespace idocp {

inline ImpulseBackwardRiccatiRecursionFactorizer::
ImpulseBackwardRiccatiRecursionFactorizer(const Robot& robot) 
  : has_floating_base_(robot.hasFloatingBase()),
    dimv_(robot.dimv()),
    AtPqq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPqv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvq_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())),
    AtPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
}


inline ImpulseBackwardRiccatiRecursionFactorizer::
ImpulseBackwardRiccatiRecursionFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    AtPqq_(),
    AtPqv_(),
    AtPvq_(),
    AtPvv_() {
}


inline ImpulseBackwardRiccatiRecursionFactorizer::
~ImpulseBackwardRiccatiRecursionFactorizer() {
}


inline void ImpulseBackwardRiccatiRecursionFactorizer::factorizeKKTMatrix(
    const SplitRiccatiFactorization& riccati_next, 
    ImpulseSplitKKTMatrix& kkt_matrix) {
  if (has_floating_base_) {
    AtPqq_.template topRows<6>().noalias() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.Pqq.template topRows<6>();
    AtPqq_.bottomRows(dimv_-6) = riccati_next.Pqq.bottomRows(dimv_-6);
    AtPqv_.template topRows<6>().noalias() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.Pqv.template topRows<6>();
    AtPqv_.bottomRows(dimv_-6) = riccati_next.Pqv.bottomRows(dimv_-6);
  }
  else {
    AtPqq_.noalias() = riccati_next.Pqq;
    AtPqv_.noalias() = riccati_next.Pqv;
  }
  AtPqq_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pqv.transpose();
  AtPqv_.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.Pvv;
  AtPvq_.noalias()  = kkt_matrix.Fvv().transpose() * riccati_next.Pqv.transpose();
  AtPvv_.noalias()  = kkt_matrix.Fvv().transpose() * riccati_next.Pvv;
  // Factorize F
  if (has_floating_base_) {
    kkt_matrix.Qqq().template leftCols<6>().noalias() 
        += AtPqq_.template leftCols<6>() 
            * kkt_matrix.Fqq().template topLeftCorner<6, 6>();
    kkt_matrix.Qqq().rightCols(dimv_-6).noalias() += AtPqq_.rightCols(dimv_-6);
  }
  else {
    kkt_matrix.Qqq().noalias() += AtPqq_;
  }
  kkt_matrix.Qqq().noalias() += AtPqv_ * kkt_matrix.Fvq();
  kkt_matrix.Qqv().noalias() += AtPqv_ * kkt_matrix.Fvv();
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qvv().noalias() += AtPvv_ * kkt_matrix.Fvv();
}


inline void ImpulseBackwardRiccatiRecursionFactorizer::
factorizeRiccatiFactorization(const SplitRiccatiFactorization& riccati_next, 
                              const ImpulseSplitKKTMatrix& kkt_matrix, 
                              const ImpulseSplitKKTResidual& kkt_residual, 
                              SplitRiccatiFactorization& riccati) {
  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pvv = kkt_matrix.Qvv();
  riccati.Pvq = riccati.Pqv.transpose();
  // preserve the symmetry
  riccati.Pqq = 0.5 * (riccati.Pqq + riccati.Pqq.transpose()).eval();
  riccati.Pvv = 0.5 * (riccati.Pvv + riccati.Pvv.transpose()).eval();
  if (has_floating_base_) {
    riccati.sq.template head<6>().noalias() 
        = kkt_matrix.Fqq().template topLeftCorner<6, 6>().transpose() 
            * riccati_next.sq.template head<6>();
    riccati.sq.tail(dimv_-6) = riccati_next.sq.tail(dimv_-6);
  }
  else {
    riccati.sq.noalias() = riccati_next.sq;
  }
  riccati.sq.noalias() += kkt_matrix.Fvq().transpose() * riccati_next.sv;
  riccati.sv.noalias()  = kkt_matrix.Fvv().transpose() * riccati_next.sv;
  riccati.sq.noalias() -= AtPqq_ * kkt_residual.Fq();
  riccati.sq.noalias() -= AtPqv_ * kkt_residual.Fv();
  riccati.sv.noalias() -= AtPvq_ * kkt_residual.Fq();
  riccati.sv.noalias() -= AtPvv_ * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();
  riccati.sv.noalias() -= kkt_residual.lv();
}

} // namespace idocp

#endif // IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HXX_