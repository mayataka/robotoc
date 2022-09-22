#include "robotoc/dynamics/impulse_dynamics.hpp"

#include <cassert>


namespace robotoc {

ImpulseDynamics::ImpulseDynamics(const Robot& robot) 
  : data_(robot) {
}


ImpulseDynamics::ImpulseDynamics() 
  : data_() {
}


void ImpulseDynamics::evalImpulseDynamics(Robot& robot, 
                                          const ImpulseStatus& impulse_status, 
                                          const SplitSolution& s) {
  data_.setContactDimension(impulse_status.dimf());
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_.ID_full());
  robot.computeImpulseVelocityResidual(impulse_status, data_.C());
}


void ImpulseDynamics::linearizeImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status,  
    const SplitSolution& s, SplitKKTResidual& kkt_residual) {
  evalImpulseDynamics(robot, impulse_status, s);
  robot.RNEAImpulseDerivatives(s.q, s.dv, data_.dIDdq(), data_.dIDddv);
  robot.computeImpulseVelocityDerivatives(impulse_status, data_.dCdq(), 
                                          data_.dCdv());
  // augment inverse impulse dynamics constraint
  kkt_residual.lq().noalias() += data_.dIDdq().transpose() * s.beta;
  kkt_residual.ldv.noalias()  += data_.dIDddv.transpose() * s.beta;
  kkt_residual.lf().noalias() -= data_.dCdv() * s.beta;
  // augment impulse velocity constraint
  kkt_residual.lq().noalias() += data_.dCdq().transpose() * s.mu_stack();
  kkt_residual.lv().noalias() += data_.dCdv().transpose() * s.mu_stack();
  kkt_residual.ldv.noalias()  += data_.dCdv().transpose() * s.mu_stack();
}


void ImpulseDynamics::condenseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  robot.computeMJtJinv(data_.dIDddv, data_.dCdv(), data_.MJtJinv());
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimf();
  data_.MJtJinv_dIDCdqv().leftCols(dimv).noalias() 
      = data_.MJtJinv() * data_.dIDCdqv().leftCols(dimv);
  data_.MJtJinv_dIDCdqv().topRightCorner(dimv, dimv).noalias() 
      = data_.MJtJinv().topRightCorner(dimv, dimf) * data_.dCdv();
  data_.MJtJinv_dIDCdqv().bottomRightCorner(dimf, dimv).noalias() 
      = data_.MJtJinv().bottomRightCorner(dimf, dimf) * data_.dCdv();
  data_.MJtJinv_IDC().noalias() = data_.MJtJinv() * data_.IDC();

  data_.Qdvfqv().topRows(dimv).noalias() 
      = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() 
          * data_.MJtJinv_dIDCdqv().topRows(dimv);
  data_.Qdvfqv().bottomRows(dimf).noalias() 
      = - kkt_matrix.Qff() * data_.MJtJinv_dIDCdqv().bottomRows(dimf);
  data_.Qdvfqv().bottomLeftCorner(dimf, dimv).noalias() 
      -= kkt_matrix.Qqf().transpose();
  data_.ldv() = kkt_residual.ldv;
  data_.lf()  = - kkt_residual.lf();
  data_.ldv().noalias() 
      -= kkt_matrix.Qdvdv.diagonal().asDiagonal() 
          * data_.MJtJinv_IDC().head(dimv);
  data_.lf().noalias() -= kkt_matrix.Qff() * data_.MJtJinv_IDC().tail(dimf);

  kkt_matrix.Qxx.noalias() 
      -= data_.MJtJinv_dIDCdqv().transpose() * data_.Qdvfqv();
  kkt_matrix.Qxx.topRows(dimv).noalias() 
      += kkt_matrix.Qqf() * data_.MJtJinv_dIDCdqv().bottomRows(dimf);
  kkt_residual.lx.noalias() 
      -= data_.MJtJinv_dIDCdqv().transpose() * data_.ldvf();
  kkt_residual.lq().noalias()
      += kkt_matrix.Qqf() * data_.MJtJinv_IDC().tail(dimf);

  kkt_matrix.Fvq() = - data_.MJtJinv_dIDCdqv().topLeftCorner(dimv, dimv);
  kkt_matrix.Fvv() = Eigen::MatrixXd::Identity(dimv, dimv) 
                    - data_.MJtJinv_dIDCdqv().topRightCorner(dimv, dimv);
  kkt_residual.Fv().noalias() -= data_.MJtJinv_IDC().head(dimv);
}


void ImpulseDynamics::expandPrimal(SplitDirection& d) const {
  d.ddvf().noalias()  = - data_.MJtJinv_dIDCdqv() * d.dx;
  d.ddvf().noalias() -= data_.MJtJinv_IDC();
  d.df().array()     *= -1;
}


void ImpulseDynamics::expandDual(const SplitDirection& d_next, SplitDirection& d) {
  data_.ldvf().noalias() += data_.Qdvfqv() * d.dx;
  data_.ldv().noalias()  += d_next.dgmm();
  d.dbetamu().noalias()   = - data_.MJtJinv() * data_.ldvf();
}

} // namespace robotoc 