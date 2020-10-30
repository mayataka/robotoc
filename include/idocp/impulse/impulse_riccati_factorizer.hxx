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
    AtPvv_(Eigen::MatrixXd::Zero(robot.dimv(), robot.dimv())) {
}


inline ImpulseRiccatiFactorizer::ImpulseRiccatiFactorizer() 
  : has_floating_base_(false),
    dimv_(0),
    AtPqq_(),
    AtPqv_(),
    AtPvq_(),
    AtPvv_() {
}


inline ImpulseRiccatiFactorizer::~ImpulseRiccatiFactorizer() {
}


inline void ImpulseRiccatiFactorizer::factorize(
    const ImpulseKKTMatrix& kkt_matrix, const ImpulseKKTResidual& kkt_residual,
    const RiccatiSolution& riccati_next, RiccatiSolution& riccati) {
  AtPqq_.noalias() = kkt_matrix.Fqq.transpose() * riccati_next.Pqq;
  AtPqq_.noalias() += kkt_matrix.Fvq.transpose() * riccati_next.Pvq;
  AtPqv_.noalias() = kkt_matrix.Fqq.transpose() * riccati_next.Pqv;
  AtPqv_.noalias() += kkt_matrix.Fvq.transpose() * riccati_next.Pvv;
  AtPvq_.noalias() = kkt_matrix.Fvv.transpose() * riccati_next.Pvq;
  AtPvv_.noalias() = kkt_matrix.Fvv.transpose() * riccati_next.Pvv;

  riccati.Pqq = kkt_matrix.Qqq();
  riccati.Pqq.noalias() += AtPqq_ * kkt_matrix.Fqq;
  riccati.Pqq.noalias() += AtPqv_ * kkt_matrix.Fvq;

  riccati.Pqv = kkt_matrix.Qqv();
  riccati.Pqv.noalias() += AtPqv_ * kkt_matrix.Fvv;

  riccati.Pvq = riccati.Pqv.transpose();

  riccati.Pvv = kkt_matrix.Qvv();
  riccati.Pvv.noalias() += AtPvv_ * kkt_matrix.Fvv;

  riccati.sq.noalias() = kkt_matrix.Fqq.transpose() * riccati_next.sq;
  riccati.sq.noalias() += kkt_matrix.Fvq.transpose() * riccati_next.sv;
  riccati.sq.noalias() -= AtPqq_ * kkt_residual.Fq();
  riccati.sq.noalias() -= AtPqv_ * kkt_residual.Fv();
  riccati.sq.noalias() -= kkt_residual.lq();

  riccati.sv.noalias() = kkt_matrix.Fvv.transpose() * riccati_next.sv;
  riccati.sv.noalias() -= AtPvq_ * kkt_residual.Fq();
  riccati.sv.noalias() -= AtPvv_ * kkt_residual.Fv();
  riccati.sv.noalias() -= kkt_residual.lv();
}

} // namespace idocp

#endif // IDOCP_IMPULSE_RICCATI_FACTORIZER_HXX_