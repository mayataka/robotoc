#ifndef IDOCP_SPLIT_BACKWARD_CORRECTION_HXX_
#define IDOCP_SPLIT_BACKWARD_CORRECTION_HXX_

#include "idocp/ocp/split_backward_correction.hpp"

#include <cassert>

namespace idocp {

inline SplitBackwardCorrection::SplitBackwardCorrection(const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimu_(robot.dimu()),
    dimKKT_(4*robot.dimv()+robot.dimu()),
    is_impulse_condition_valid_(false),
    kkt_mat_inverter_(robot),
    data_(robot),
    x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())) {
}


inline SplitBackwardCorrection::SplitBackwardCorrection() 
  : dimv_(0),
    dimx_(0),
    dimu_(0),
    dimKKT_(0),
    is_impulse_condition_valid_(false),
    kkt_mat_inverter_(),
    data_(),
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
    const SplitSolution& s, SplitSolution& s_new) {
  assert(aux_mat_next.rows() == dimx_);
  assert(aux_mat_next.cols() == dimx_);
  kkt_matrix.Qxx().noalias() += aux_mat_next;
  coarseUpdate(robot, dtau, kkt_matrix, kkt_residual, s, s_new);
}


inline void SplitBackwardCorrection::coarseUpdate(
    const Robot& robot, const double dtau, SplitKKTMatrix& kkt_matrix, 
    const SplitKKTResidual& kkt_residual, const SplitSolution& s, 
    SplitSolution& s_new) {
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qux() = kkt_matrix.Qxu().transpose();
  data_.setImpulseStatus(0);
  dimKKT_ = 2*dimx_ + dimu_;
  is_impulse_condition_valid_ = false;
  kkt_mat_inverter_.invert(dtau, kkt_matrix.Jac(), kkt_matrix.Qss(), 
                           data_.KKT_mat_inv());
  data_.splitDirection().noalias() 
      = data_.KKT_mat_inv() * kkt_residual.splitKKTResidual();
  s_new.lmd = s.lmd - data_.dlmd();
  s_new.gmm = s.gmm - data_.dgmm();
  s_new.u   = s.u   - data_.du(); 
  robot.integrateConfiguration(s.q, data_.dq(), -1, s_new.q);
  s_new.v   = s.v   - data_.dv();
}


template <typename MatrixType>
inline void SplitBackwardCorrection::coarseUpdate(
    const Robot& robot, const double dtau, 
    const Eigen::MatrixBase<MatrixType>& aux_mat_next,
    SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
    const SplitSolution& s, const ImpulseSplitSolution& s_next, 
    SplitSolution& s_new) {
  assert(aux_mat_next.rows() == dimx_);
  assert(aux_mat_next.cols() == dimx_);
  kkt_matrix.Qxx().noalias() += aux_mat_next;
  kkt_matrix.Qvq() = kkt_matrix.Qqv().transpose();
  kkt_matrix.Qux() = kkt_matrix.Qxu().transpose();
  data_.setImpulseStatus(s_next.dimf());
  dimKKT_ = 2*dimx_ + dimu_ + s_next.dimf();
  is_impulse_condition_valid_ = true;
  kkt_mat_inverter_.invert(dtau, kkt_matrix.Jac(), kkt_matrix.Pq(), 
                           kkt_matrix.Qss(), data_.KKT_mat_inv());
  data_.splitDirection().noalias() 
      = data_.KKT_mat_inv() * kkt_residual.splitKKTResidual();
  s_new.lmd        = s.lmd - data_.dlmd();
  s_new.gmm        = s.gmm - data_.dgmm();
  data_.xi_stack() = s_next.xi_stack() - data_.dxi();
  s_new.u          = s.u - data_.du(); 
  robot.integrateConfiguration(s.q, data_.dq(), -1, s_new.q);
  s_new.v          = s.v - data_.dv();
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitBackwardCorrection::auxMat() const {
  return data_.auxMat();
}


inline void SplitBackwardCorrection::backwardCorrectionSerial(
    const SplitSolution& s_next, const SplitSolution& s_new_next,
    SplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(dimv_) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() 
      = data_.KKT_mat_inv().block(0, dimKKT_-dimx_, dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(dimv_);
  s_new.gmm.noalias() -= dx_.tail(dimv_);
}


inline void SplitBackwardCorrection::backwardCorrectionParallel(
    const Robot& robot, SplitSolution& s_new) {
  data_.splitDirection().tail(dimKKT_-dimx_).noalias()
      = data_.KKT_mat_inv().block(dimx_, dimKKT_-dimx_, dimKKT_-dimx_, dimx_) 
          * x_res_;
  if (is_impulse_condition_valid_) {
    data_.xi_stack().noalias() -= data_.dxi();
  }
  s_new.u.noalias() -= data_.du();
  robot.integrateConfiguration(data_.dq(), -1, s_new.q);
  s_new.v.noalias() -= data_.dv();
}


inline void SplitBackwardCorrection::forwardCorrectionSerial(
    const Robot& robot, const SplitSolution& s_prev, 
    const SplitSolution& s_new_prev, SplitSolution& s_new) {
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, x_res_.head(dimv_));
  x_res_.tail(dimv_) = s_new_prev.v - s_prev.v;
  dx_.noalias() 
      = data_.KKT_mat_inv().block(dimKKT_-dimx_, 0, dimx_, dimx_) * x_res_;
  robot.integrateConfiguration(dx_.head(dimv_), -1, s_new.q);
  s_new.v.noalias() -= dx_.tail(dimv_);
}


inline void SplitBackwardCorrection::forwardCorrectionParallel(
    SplitSolution& s_new) {
  data_.splitDirection().head(dimKKT_-dimx_).noalias()
      = data_.KKT_mat_inv().topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= data_.dlmd();
  s_new.gmm.noalias() -= data_.dgmm();
  if (is_impulse_condition_valid_) {
    data_.xi_stack().noalias() -= data_.dxi();
  }
  s_new.u.noalias()   -= data_.du();
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


inline void SplitBackwardCorrection::computeDirection(
    const Robot& robot, const SplitSolution& s, 
    const ImpulseSplitSolution& s_next, const SplitSolution& s_new,
    SplitDirection& d, ImpulseSplitDirection& d_next) const {
  computeDirection(robot, s, s_new, d);
  d_next.dxi() = data_.xi_stack() - s_next.xi_stack();
}

} // namespace idocp

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_HXX_ 