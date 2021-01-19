#ifndef IDOCP_SPLIT_BACKWARD_CORRECTION_HXX_
#define IDOCP_SPLIT_BACKWARD_CORRECTION_HXX_

#include "idocp/ocp/split_backward_correction.hpp"

#include <cassert>

namespace idocp {

inline SplitBackwardCorrection::SplitBackwardCorrection(const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimKKT_(4*robot.dimv()+robot.dimu()),
    kkt_mat_inverter_(robot),
    KKT_mat_inv_(Eigen::MatrixXd::Zero(4*robot.dimv()+robot.dimu(), 
                                       4*robot.dimv()+robot.dimu())),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


inline SplitBackwardCorrection::SplitBackwardCorrection() 
  : dimv_(),
    dimx_(),
    dimKKT_(),
    kkt_mat_inverter_(),
    KKT_mat_inv_(),
    x_res_(),
    dx_() {
}


inline SplitBackwardCorrection::~SplitBackwardCorrection() {
}


template <typename MatrixType>
inline void SplitBackwardCorrection::coarseUpdate(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
    SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const SplitSolution& s, SplitDirection& d, SplitSolution& s_new) {
  assert(aux_mat_next.rows() == dimx_);
  assert(aux_mat_next.cols() == dimx_);
  kkt_matrix.Qxx().noalias() += aux_mat_next;
  coarseUpdate(robot, dtau, kkt_matrix, kkt_residual, s, d, s_new);
}


inline void SplitBackwardCorrection::coarseUpdate(
    const Robot& robot, const double dtau, SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, const SplitSolution& s, 
    SplitDirection& d, SplitSolution& s_new) {
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qux() = kkt_matrix.Qxu().transpose();
  kkt_mat_inverter_.invert(dtau, kkt_matrix.Jac(), kkt_matrix.Qss(), KKT_mat_inv_);
  d.splitDirection().noalias() = KKT_mat_inv_ * kkt_residual.splitKKTResidual();
  s_new.lmd = s.lmd - d.dlmd();
  s_new.gmm = s.gmm - d.dgmm();
  s_new.u   = s.u - d.du(); 
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new.q);
  s_new.v   = s.v - d.dv();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitBackwardCorrection::auxMat() const {
  return KKT_mat_inv_.topLeftCorner(dimx_, dimx_);
}


inline void SplitBackwardCorrection::backwardCorrectionSerial(
    const SplitSolution& s_next, const SplitSolution& s_new_next,
    SplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(dimv_) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() = KKT_mat_inv_.block(0, dimKKT_-dimx_, dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(dimv_);
  s_new.gmm.noalias() -= dx_.tail(dimv_);
}


inline void SplitBackwardCorrection::backwardCorrectionParallel(
    const Robot& robot, SplitDirection& d, SplitSolution& s_new) const {
  d.splitDirection().tail(dimKKT_-dimx_).noalias()
      = KKT_mat_inv_.block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) * x_res_;
  s_new.u.noalias() -= d.du();
  robot.integrateConfiguration(d.dq(), -1, s_new.q);
  s_new.v.noalias() -= d.dv();
}


inline void SplitBackwardCorrection::forwardCorrectionSerial(
    const Robot& robot, const SplitSolution& s_prev, 
    const SplitSolution& s_new_prev, SplitSolution& s_new) {
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, x_res_.head(dimv_));
  x_res_.tail(dimv_) = s_new_prev.v - s_prev.v;
  dx_.noalias() = KKT_mat_inv_.block(dimKKT_-dimx_, 0, dimx_, dimx_) * x_res_;
  robot.integrateConfiguration(dx_.head(dimv_), -1, s_new.q);
  s_new.v.noalias() -= dx_.tail(dimv_);
}


inline void SplitBackwardCorrection::forwardCorrectionParallel(
    SplitDirection& d, SplitSolution& s_new) const {
  d.splitDirection().head(dimKKT_-dimx_).noalias()
      = KKT_mat_inv_.topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= d.dlmd();
  s_new.gmm.noalias() -= d.dgmm();
  s_new.u.noalias()   -= d.du();
}


inline void SplitBackwardCorrection::computeDirection(
    const Robot& robot, const SplitSolution& s, const SplitSolution& s_new, 
    SplitDirection& d) {
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.du()   = s_new.u - s.u;
  robot.subtractConfiguration(s_new.q, s.q, d.dq());
  d.dv()   = s_new.v - s.v;
}

} // namespace idocp

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_HXX_ 