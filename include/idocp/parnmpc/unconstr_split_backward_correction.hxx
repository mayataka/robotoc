#ifndef IDOCP_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HXX_
#define IDOCP_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HXX_

#include "idocp/parnmpc/unconstr_split_backward_correction.hxx"

#include <cassert>

namespace idocp {

inline UnconstrSplitBackwardCorrection::UnconstrSplitBackwardCorrection(
    const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimkkt_(5*robot.dimv()),
    kkt_mat_inverter_(robot),
    H_(Eigen::MatrixXd::Zero(3*robot.dimv(), 3*robot.dimv())),
    kkt_mat_inv_(Eigen::MatrixXd::Zero(5*robot.dimv(), 5*robot.dimv())),
    kkt_res_(Eigen::VectorXd::Zero(5*robot.dimv())),
    d_(Eigen::VectorXd::Zero(5*robot.dimv())),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


inline UnconstrSplitBackwardCorrection::UnconstrSplitBackwardCorrection() 
  : dimv_(),
    dimx_(),
    dimkkt_(),
    kkt_mat_inverter_(),
    H_(),
    kkt_mat_inv_(),
    kkt_res_(),
    d_(),
    x_res_(),
    dx_() {
}


inline UnconstrSplitBackwardCorrection::~UnconstrSplitBackwardCorrection() {
}


template <typename MatrixType>
inline void UnconstrSplitBackwardCorrection::coarseUpdate(
    const Eigen::MatrixBase<MatrixType>& aux_mat_next, const double dt,
    const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const SplitSolution& s, SplitSolution& s_new) {
  assert(aux_mat_next.rows() == dimx_);
  assert(aux_mat_next.cols() == dimx_);
  H_.topLeftCorner(dimv_, dimv_)     = kkt_matrix.Qaa;
  H_.bottomLeftCorner(dimx_, dimv_)  = kkt_matrix.Qxu; // This is actually Qxa
  H_.bottomRightCorner(dimx_, dimx_) = kkt_matrix.Qxx.transpose();
  H_.bottomRightCorner(dimx_, dimx_).noalias() += aux_mat_next;
  kkt_mat_inverter_.invert(dt, H_, kkt_mat_inv_);
  kkt_res_.head(dimx_)           = kkt_residual.Fx;
  kkt_res_.segment(dimx_, dimv_) = kkt_residual.la;
  kkt_res_.tail(dimx_)           = kkt_residual.lx;
  d_.noalias() = kkt_mat_inv_ * kkt_res_;

  s_new.lmd = s.lmd - d_.head(dimv_);
  s_new.gmm = s.gmm - d_.segment(dimv_, dimv_);
  s_new.a   = s.a   - d_.segment(2*dimv_, dimv_); 
  s_new.q   = s.q   - d_.segment(3*dimv_, dimv_);
  s_new.v   = s.v   - d_.tail(dimv_);
}


inline void UnconstrSplitBackwardCorrection::coarseUpdate(
    const double dt, const SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, const SplitSolution& s, 
    SplitSolution& s_new) {
  H_.topLeftCorner(dimv_, dimv_)     = kkt_matrix.Qaa;
  H_.bottomLeftCorner(dimx_, dimv_)  = kkt_matrix.Qxu; // This is actually Qxa
  H_.bottomRightCorner(dimx_, dimx_) = kkt_matrix.Qxx.transpose();
  kkt_mat_inverter_.invert(dt, H_, kkt_mat_inv_);
  kkt_res_.head(dimx_)           = kkt_residual.Fx;
  kkt_res_.segment(dimx_, dimv_) = kkt_residual.la;
  kkt_res_.tail(dimx_)           = kkt_residual.lx;
  d_.noalias() = kkt_mat_inv_ * kkt_res_;

  s_new.lmd = s.lmd - d_.head(dimv_);
  s_new.gmm = s.gmm - d_.segment(dimv_, dimv_);
  s_new.a   = s.a   - d_.segment(2*dimv_, dimv_); 
  s_new.q   = s.q   - d_.segment(3*dimv_, dimv_);
  s_new.v   = s.v   - d_.tail(dimv_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
UnconstrSplitBackwardCorrection::auxMat() const {
  return kkt_mat_inv_.topLeftCorner(dimx_, dimx_);
}


inline void UnconstrSplitBackwardCorrection::backwardCorrectionSerial(
    const SplitSolution& s_next, const SplitSolution& s_new_next,
    SplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(dimv_) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() = kkt_mat_inv_.block(0, dimkkt_-dimx_, dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(dimv_);
  s_new.gmm.noalias() -= dx_.tail(dimv_);
}


inline void UnconstrSplitBackwardCorrection::backwardCorrectionParallel(
    SplitSolution& s_new) {
  d_.tail(dimkkt_-dimx_).noalias()
      = kkt_mat_inv_.block(dimx_, dimkkt_-dimx_, dimkkt_-dimx_, dimx_) * x_res_;
  s_new.a.noalias() -= d_.segment(2*dimv_, dimv_);
  s_new.q.noalias() -= d_.segment(3*dimv_, dimv_);
  s_new.v.noalias() -= d_.tail(dimv_);
}


inline void UnconstrSplitBackwardCorrection::forwardCorrectionSerial(
    const SplitSolution& s_prev, const SplitSolution& s_new_prev, 
    SplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_prev.q - s_prev.q;
  x_res_.tail(dimv_) = s_new_prev.v - s_prev.v;
  dx_.noalias() = kkt_mat_inv_.block(dimkkt_-dimx_, 0, dimx_, dimx_) * x_res_;
  s_new.q.noalias() -= dx_.head(dimv_);
  s_new.v.noalias() -= dx_.tail(dimv_);
}


inline void UnconstrSplitBackwardCorrection::forwardCorrectionParallel(
    SplitSolution& s_new) {
  d_.head(dimkkt_-dimx_).noalias()
      = kkt_mat_inv_.topLeftCorner(dimkkt_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= d_.head(dimv_);
  s_new.gmm.noalias() -= d_.segment(dimv_, dimv_);
  s_new.a.noalias()   -= d_.segment(2*dimv_, dimv_);
}


inline void UnconstrSplitBackwardCorrection::computeDirection(
    const SplitSolution& s, const SplitSolution& s_new, 
    SplitDirection& d) {
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.da()   = s_new.a - s.a;
  d.dq()   = s_new.q - s.q;
  d.dv()   = s_new.v - s.v;
}

} // namespace idocp

#endif // IDOCP_UNCONSTR_SPLIT_BACKWARD_CORRECTION_HXX_ 