#ifndef IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HXX_
#define IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HXX_

#include "idocp/impulse/impulse_split_backward_correction.hpp"

#include <cassert>

namespace idocp {

inline ImpulseSplitBackwardCorrection::ImpulseSplitBackwardCorrection(
      const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimf_(0),
    dimKKT_(5*robot.dimv()),
    kkt_mat_inverter_(robot.dimv(), robot.max_dimf()),
    KKT_mat_inv_(Eigen::MatrixXd::Zero(4*robot.dimv()+2*robot.max_dimf(), 
                                       4*robot.dimv()+2*robot.max_dimf())),
    split_direction_full_(
        Eigen::VectorXd::Zero(4*robot.dimv()+2*robot.max_dimf())),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


inline ImpulseSplitBackwardCorrection::ImpulseSplitBackwardCorrection() 
  : dimv_(0),
    dimx_(0),
    dimf_(0),
    dimKKT_(0),
    kkt_mat_inverter_(),
    KKT_mat_inv_(),
    split_direction_full_(),
    x_res_(),
    dx_() {
}


inline ImpulseSplitBackwardCorrection::~ImpulseSplitBackwardCorrection() {
}


template <typename MatrixType>
inline void ImpulseSplitBackwardCorrection::coarseUpdate( 
    const Eigen::MatrixBase<MatrixType>& aux_mat_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const ImpulseSplitSolution& s, ImpulseSplitSolution& s_new) {
  assert(aux_mat_next.rows() == dimx_);
  assert(aux_mat_next.cols() == dimx_);
  kkt_matrix.Qxx().noalias() += aux_mat_next;
  dimf_ = s.dimf();
  dimKKT_ = kkt_matrix.dimKKT();
  kkt_mat_inverter_.invert(kkt_matrix.Jac(), kkt_matrix.Qss(), 
                           KKT_mat_inv_.topLeftCorner(dimKKT_, dimKKT_));
  split_direction().noalias() = KKT_mat_inv_.topLeftCorner(dimKKT_, dimKKT_) 
                                  * kkt_residual.KKT_residual();
  s_new.lmd        = s.lmd         - dlmd();
  s_new.gmm        = s.gmm         - dgmm();
  s_new.mu_stack() = s.mu_stack()  - dmu(); 
  s_new.f_stack()  = s.f_stack()   - df(); 
  s_new.q          = s.q           - dq();
  s_new.v          = s.v           - dv();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitBackwardCorrection::auxMat() const {
  return KKT_mat_inv_.topLeftCorner(dimx_, dimx_);
}


inline void ImpulseSplitBackwardCorrection::backwardCorrectionSerial(
    const SplitSolution& s_next, const SplitSolution& s_new_next,
    ImpulseSplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(dimv_) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() = KKT_mat_inv_.block(0, dimKKT_-dimx_, dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(dimv_);
  s_new.gmm.noalias() -= dx_.tail(dimv_);
}


inline void ImpulseSplitBackwardCorrection::backwardCorrectionParallel(
    ImpulseSplitSolution& s_new) {
  split_direction().tail(dimKKT_-dimx_).noalias()
      = KKT_mat_inv_.block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) * x_res_;
  s_new.mu_stack().noalias() -= dmu();
  s_new.f_stack().noalias()  -= df();
  s_new.q.noalias()          -= dq();
  s_new.v.noalias()          -= dv();
}


inline void ImpulseSplitBackwardCorrection::forwardCorrectionSerial(
    const SplitSolution& s_prev, const SplitSolution& s_new_prev, 
    ImpulseSplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_prev.q - s_prev.q;
  x_res_.tail(dimv_) = s_new_prev.v - s_prev.v;
  dx_.noalias() = KKT_mat_inv_.block(dimKKT_-dimx_, 0, dimx_, dimx_) * x_res_;
  s_new.q.noalias() -= dx_.head(dimv_);
  s_new.v.noalias() -= dx_.tail(dimv_);
}


inline void ImpulseSplitBackwardCorrection::forwardCorrectionParallel(
    ImpulseSplitSolution& s_new) {
  split_direction().head(dimKKT_-dimx_).noalias()
      = KKT_mat_inv_.topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias()        -= dlmd();
  s_new.gmm.noalias()        -= dgmm();
  s_new.mu_stack().noalias() -= dmu();
  s_new.f_stack().noalias()  -= df();
}


inline void ImpulseSplitBackwardCorrection::computeDirection(
    const ImpulseSplitSolution& s, const ImpulseSplitSolution& s_new, 
    ImpulseSplitDirection& d) {
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.dmu()  = s_new.mu_stack() - s.mu_stack();
  d.df()   = s_new.f_stack() - s.f_stack();
  d.dq()   = s_new.q - s.q;
  d.dv()   = s_new.v - s.v;
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrection::split_direction() {
  return split_direction_full_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrection::dlmd() {
  return split_direction_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrection::dgmm() {
  return split_direction_full_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrection::dmu() {
  return split_direction_full_.segment(2*dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrection::df() {
  return split_direction_full_.segment(2*dimv_+dimf_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrection::dq() {
  return split_direction_full_.segment(2*dimv_+2*dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrection::dv() {
  return split_direction_full_.segment(3*dimv_+2*dimf_, dimv_);
}

} // namespace idocp

#endif // IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HXX_ 