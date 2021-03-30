#ifndef IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HXX_
#define IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HXX_

#include "idocp/impulse/impulse_split_backward_correction.hpp"

#include <cassert>

namespace idocp {

inline ImpulseSplitBackwardCorrection::ImpulseSplitBackwardCorrection(
      const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimKKT_(4*robot.dimv()),
    kkt_mat_inverter_(robot),
    data_(robot),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


inline ImpulseSplitBackwardCorrection::ImpulseSplitBackwardCorrection() 
  : dimv_(0),
    dimx_(0),
    dimKKT_(0),
    kkt_mat_inverter_(),
    data_(),
    x_res_(),
    dx_() {
}


inline ImpulseSplitBackwardCorrection::~ImpulseSplitBackwardCorrection() {
}


template <typename MatrixType>
inline void ImpulseSplitBackwardCorrection::coarseUpdate( 
    const Robot& robot, const Eigen::MatrixBase<MatrixType>& aux_mat_next, 
    ImpulseSplitKKTMatrix& kkt_matrix, 
    const ImpulseSplitKKTResidual& kkt_residual, 
    const ImpulseSplitSolution& s, ImpulseSplitSolution& s_new) {
  assert(aux_mat_next.rows() == dimx_);
  assert(aux_mat_next.cols() == dimx_);
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  data_.setImpulseStatus(s.dimf());
  dimKKT_ = 2*dimx_ + 2*s.dimf();
  kkt_matrix.Qxx().noalias() += aux_mat_next;
  kkt_mat_inverter_.invert(kkt_matrix.Jac(), kkt_matrix.Qss(), 
                           data_.KKT_mat_inv());
  data_.splitDirection().noalias() 
      = data_.KKT_mat_inv() * kkt_residual.splitKKTResidual();
  s_new.setImpulseStatus(s);
  s_new.lmd        = s.lmd         - data_.dlmd();
  s_new.gmm        = s.gmm         - data_.dgmm();
  s_new.mu_stack() = s.mu_stack()  - data_.dmu(); 
  s_new.f_stack()  = s.f_stack()   - data_.df(); 
  robot.integrateConfiguration(s.q,  data_.dq(), -1, s_new.q);
  s_new.v          = s.v           - data_.dv();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitBackwardCorrection::auxMat() const {
  return data_.auxMat();
}


inline void ImpulseSplitBackwardCorrection::backwardCorrectionSerial(
    const SplitSolution& s_next, const SplitSolution& s_new_next,
    ImpulseSplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(dimv_) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() 
      = data_.KKT_mat_inv().block(0, dimKKT_-dimx_, dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(dimv_);
  s_new.gmm.noalias() -= dx_.tail(dimv_);
}


inline void ImpulseSplitBackwardCorrection::backwardCorrectionParallel(
    const Robot& robot, ImpulseSplitSolution& s_new) {
  data_.splitDirection().tail(dimKKT_-dimx_).noalias()
      = data_.KKT_mat_inv().block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) 
          * x_res_;
  s_new.mu_stack().noalias() -= data_.dmu();
  s_new.f_stack().noalias()  -= data_.df();
  robot.integrateConfiguration(data_.dq(), -1, s_new.q);
  s_new.v.noalias()          -= data_.dv();
}


inline void ImpulseSplitBackwardCorrection::forwardCorrectionSerial(
    const Robot& robot, const SplitSolution& s_prev, 
    const SplitSolution& s_new_prev, ImpulseSplitSolution& s_new) {
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, x_res_.head(dimv_));
  x_res_.tail(dimv_) = s_new_prev.v - s_prev.v;
  dx_.noalias() 
      = data_.KKT_mat_inv().block(dimKKT_-dimx_, 0, dimx_, dimx_) * x_res_;
  robot.integrateConfiguration(dx_.head(dimv_), -1, s_new.q);
  s_new.v.noalias() -= dx_.tail(dimv_);
}


inline void ImpulseSplitBackwardCorrection::forwardCorrectionParallel(
    ImpulseSplitSolution& s_new) {
  data_.splitDirection().head(dimKKT_-dimx_).noalias()
      = data_.KKT_mat_inv().topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias()        -= data_.dlmd();
  s_new.gmm.noalias()        -= data_.dgmm();
  s_new.mu_stack().noalias() -= data_.dmu();
  s_new.f_stack().noalias()  -= data_.df();
}


inline void ImpulseSplitBackwardCorrection::computeDirection(
    const Robot& robot, const ImpulseSplitSolution& s, 
    const ImpulseSplitSolution& s_new, ImpulseSplitDirection& d) {
  d.setImpulseStatusByDimension(s.dimf());
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.dmu()  = s_new.mu_stack() - s.mu_stack();
  d.df()   = s_new.f_stack() - s.f_stack();
  robot.subtractConfiguration(s_new.q, s.q, d.dq());
  d.dv()   = s_new.v - s.v;
}

} // namespace idocp

#endif // IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_HXX_ 