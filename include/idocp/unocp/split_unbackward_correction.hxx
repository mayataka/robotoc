#ifndef IDOCP_SPLIT_UNBACKWARD_CORRECTION_HXX_
#define IDOCP_SPLIT_UNBACKWARD_CORRECTION_HXX_

#include "idocp/unocp/split_unbackward_correction.hpp"

#include <cassert>

namespace idocp {

inline SplitUnBackwardCorrection::SplitUnBackwardCorrection(const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimKKT_(5*robot.dimv()),
    KKT_mat_inv_(Eigen::MatrixXd::Zero(5*robot.dimv(), 5*robot.dimv())),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


inline SplitUnBackwardCorrection::SplitUnBackwardCorrection() 
  : dimv_(),
    dimx_(),
    dimKKT_(),
    KKT_mat_inv_(),
    x_res_(),
    dx_() {
}


inline SplitUnBackwardCorrection::~SplitUnBackwardCorrection() {
}


template <typename MatrixType>
inline void SplitUnBackwardCorrection::coarseUpdate(
    const Eigen::MatrixBase<MatrixType>& aux_mat_next, const double dtau,
    SplitUnKKTMatrix& unkkt_matrix, const SplitUnKKTResidual& unkkt_residual, 
    const SplitSolution& s, SplitDirection& d, SplitSolution& s_new) {
  assert(aux_mat_next.rows() == dimx_);
  assert(aux_mat_next.cols() == dimx_);
  unkkt_matrix.Qxx().noalias() += aux_mat_next;
  unkkt_matrix.invert(dtau , KKT_mat_inv_);
  d.split_direction.noalias() = KKT_mat_inv_ * unkkt_residual.KKT_residual;
  s_new.lmd = s.lmd - d.dlmd();
  s_new.gmm = s.gmm - d.dgmm();
  s_new.a   = s.a - d.da();
  s_new.q   = s.q - d.dq();
  s_new.v   = s.v - d.dv();
}


inline void SplitUnBackwardCorrection::coarseUpdate(
    const double dtau, SplitUnKKTMatrix& unkkt_matrix, 
    const SplitUnKKTResidual& unkkt_residual, const SplitSolution& s, 
    SplitDirection& d, SplitSolution& s_new) {
  unkkt_matrix.invert(dtau , KKT_mat_inv_);
  d.split_direction.noalias() = KKT_mat_inv_ * unkkt_residual.KKT_residual;
  s_new.lmd = s.lmd - d.dlmd();
  s_new.gmm = s.gmm - d.dgmm();
  s_new.a   = s.a - d.da();
  s_new.q   = s.q - d.dq();
  s_new.v   = s.v - d.dv();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitUnBackwardCorrection::auxMat() const {
  return KKT_mat_inv_.topLeftCorner(dimx_, dimx_);
}


inline void SplitUnBackwardCorrection::backwardCorrectionSerial(
    const SplitSolution& s_next, const SplitSolution& s_new_next,
    SplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(dimv_) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() = KKT_mat_inv_.block(0, dimKKT_-dimx_, dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(dimv_);
  s_new.gmm.noalias() -= dx_.tail(dimv_);
}


inline void SplitUnBackwardCorrection::backwardCorrectionParallel(
    SplitDirection& d, SplitSolution& s_new) const {
  d.split_direction.tail(dimKKT_-dimx_).noalias()
      = KKT_mat_inv_.block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) * x_res_;
  s_new.a.noalias() -= d.da();
  s_new.q.noalias() -= d.dq();
  s_new.v.noalias() -= d.dv();
}


inline void SplitUnBackwardCorrection::forwardCorrectionSerial(
    const SplitSolution& s_prev, const SplitSolution& s_new_prev, 
    SplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_prev.q - s_prev.q;
  x_res_.tail(dimv_) = s_new_prev.v - s_prev.v;
  dx_.noalias() = KKT_mat_inv_.block(dimKKT_-dimx_, 0, dimx_, dimx_) * x_res_;
  s_new.q.noalias() -= dx_.head(dimv_);
  s_new.v.noalias() -= dx_.tail(dimv_);
}


inline void SplitUnBackwardCorrection::forwardCorrectionParallel(
    SplitDirection& d, SplitSolution& s_new) const {
  d.split_direction.head(dimKKT_-dimx_).noalias()
      = KKT_mat_inv_.topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= d.dlmd();
  s_new.gmm.noalias() -= d.dgmm();
  s_new.a.noalias()   -= d.da();
}


inline void SplitUnBackwardCorrection::computeDirection(
    const SplitSolution& s, const SplitSolution& s_new, 
    SplitDirection& d) {
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.da()   = s_new.a - s.a;
  d.dq()   = s_new.q - s.q;
  d.dv()   = s_new.v - s.v;
}

} // namespace idocp

#endif // IDOCP_SPLIT_UNBACKWARD_CORRECTION_HXX_ 