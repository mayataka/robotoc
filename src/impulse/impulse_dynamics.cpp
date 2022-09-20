#include "robotoc/impulse/impulse_dynamics.hpp"

#include <cassert>


namespace robotoc {

ImpulseDynamics::ImpulseDynamics(const Robot& robot) 
  : data_(robot) {
}


ImpulseDynamics::ImpulseDynamics() 
  : data_() {
}


ImpulseDynamics::~ImpulseDynamics() {
}


void ImpulseDynamics::evalImpulseDynamics(Robot& robot, 
                                          const ImpulseStatus& impulse_status, 
                                          const ImpulseSplitSolution& s) {
  data_.setImpulseStatus(impulse_status);
  robot.setImpulseForces(impulse_status, s.f);
  robot.RNEAImpulse(s.q, s.dv, data_.ImD());
  robot.computeImpulseVelocityResidual(impulse_status, data_.C());
}


void ImpulseDynamics::linearizeImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status,  
    const ImpulseSplitSolution& s, SplitKKTResidual& kkt_residual) {
  evalImpulseDynamics(robot, impulse_status, s);
  robot.RNEAImpulseDerivatives(s.q, s.dv, data_.dImDdq(), data_.dImDddv);
  robot.computeImpulseVelocityDerivatives(impulse_status, data_.dCdq(), 
                                          data_.dCdv());
  // augment inverse impulse dynamics constraint
  kkt_residual.lq().noalias() += data_.dImDdq().transpose() * s.beta;
  kkt_residual.ldv.noalias()  += data_.dImDddv.transpose() * s.beta;
  kkt_residual.lf().noalias() -= data_.dCdv() * s.beta;
  // augment impulse velocity constraint
  kkt_residual.lq().noalias() += data_.dCdq().transpose() * s.mu_stack();
  kkt_residual.lv().noalias() += data_.dCdv().transpose() * s.mu_stack();
  kkt_residual.ldv.noalias()  += data_.dCdv().transpose() * s.mu_stack();
}


void ImpulseDynamics::condenseImpulseDynamics(
    Robot& robot, const ImpulseStatus& impulse_status, 
    SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual) {
  robot.computeMJtJinv(data_.dImDddv, data_.dCdv(), data_.MJtJinv());
  const int dimv = robot.dimv();
  const int dimf = impulse_status.dimi();
  data_.MJtJinv_dImDCdqv().leftCols(dimv).noalias() 
      = data_.MJtJinv() * data_.dImDCdq();
  data_.MJtJinv_dImDCdqv().topRightCorner(dimv, dimv).noalias() 
      = data_.MJtJinv().topRightCorner(dimv, dimf) * data_.dCdv();
  data_.MJtJinv_dImDCdqv().bottomRightCorner(dimf, dimv).noalias() 
      = data_.MJtJinv().bottomRightCorner(dimf, dimf) * data_.dCdv();
  data_.MJtJinv_ImDC().noalias() = data_.MJtJinv() * data_.ImDC();

  data_.Qdvfqv().topRows(dimv).noalias() 
      = (- kkt_matrix.Qdvdv.diagonal()).asDiagonal() 
          * data_.MJtJinv_dImDCdqv().topRows(dimv);
  data_.Qdvfqv().bottomRows(dimf).noalias() 
      = - kkt_matrix.Qff() * data_.MJtJinv_dImDCdqv().bottomRows(dimf);
  data_.Qdvfqv().bottomLeftCorner(dimf, dimv).noalias() 
      -= kkt_matrix.Qqf().transpose();
  data_.ldv() = kkt_residual.ldv;
  data_.lf()  = - kkt_residual.lf();
  data_.ldv().noalias() 
      -= kkt_matrix.Qdvdv.diagonal().asDiagonal() 
          * data_.MJtJinv_ImDC().head(dimv);
  data_.lf().noalias() -= kkt_matrix.Qff() * data_.MJtJinv_ImDC().tail(dimf);

  kkt_matrix.Qxx.noalias() 
      -= data_.MJtJinv_dImDCdqv().transpose() * data_.Qdvfqv();
  kkt_matrix.Qxx.topRows(dimv).noalias() 
      += kkt_matrix.Qqf() * data_.MJtJinv_dImDCdqv().bottomRows(dimf);
  kkt_residual.lx.noalias() 
      -= data_.MJtJinv_dImDCdqv().transpose() * data_.ldvf();
  kkt_residual.lq().noalias()
      += kkt_matrix.Qqf() * data_.MJtJinv_ImDC().tail(dimf);

  kkt_matrix.Fvq() = - data_.MJtJinv_dImDCdqv().topLeftCorner(dimv, dimv);
  kkt_matrix.Fvv() = Eigen::MatrixXd::Identity(dimv, dimv) 
                    - data_.MJtJinv_dImDCdqv().topRightCorner(dimv, dimv);
  kkt_residual.Fv().noalias() -= data_.MJtJinv_ImDC().head(dimv);
}

} // namespace robotoc 