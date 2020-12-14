#ifndef IDOCP_BACKWARD_CORRECTION_HXX_
#define IDOCP_BACKWARD_CORRECTION_HXX_

#include "idocp/ocp/backward_correction.hpp"

#include <cassert>


namespace idocp {

inline BackwardCorrection::BackwardCorrection(const Robot& robot) 
  : x_res_(Eigen::VectorXd::Zero(2*robot.dimv())),
    dx_(Eigen::VectorXd::Zero(2*robot.dimv())),
    kkt_matrix_inverse_(Eigen::MatrixXd::Zero(4*robot.dimv()+robot.dimu(), 
                                              4*robot.dimv()+robot.dimu())),
    dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimKKT_(4*robot.dimv()+robot.dimu()) {
}


inline BackwardCorrection::BackwardCorrection() 
  : x_res_(),
    dx_(),
    kkt_matrix_inverse_(),
    dimv_(0),
    dimx_(0),
    dimKKT_() {
}


inline BackwardCorrection::~BackwardCorrection() {
}


inline void BackwardCorrection::coarseUpdate(const Robot& robot,
                                             const SplitSolution& s, 
                                             SplitDirection& d, 
                                             SplitKKTMatrix& kkt_matrix, 
                                             const SplitKKTResidual& kkt_residual,
                                             SplitSolution& s_new_coarse) {
  kkt_matrix.invert(kkt_matrix_inverse_);
  d.split_direction = kkt_matrix_inverse_ * kkt_residual.KKT_residual;
  s_new_coarse.lmd = s.lmd - d.dlmd();
  s_new_coarse.gmm = s.gmm - d.dgmm();
  robot.integrateConfiguration(s.q, d.dq(), -1, s_new_coarse.q);
  s_new_coarse.v = s.v - d.dv();
  s_new_coarse.u = s.u - d.du();
}


template <typename SplitSolutionType>
inline void BackwardCorrection::backwardCorrectionSerial(
    const SplitSolutionType& s_next, const SplitSolutionType& s_new_next, 
    SplitSolution& s_new) {
  x_res_.head(dimv_) = s_new_next.lmd - s_next.lmd;
  x_res_.tail(dimv_) = s_new_next.gmm - s_next.gmm;
  dx_.noalias() = kkt_matrix_inverse_.topRightCorner(dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= dx_.head(dimv_);
  s_new.gmm.noalias() -= dx_.tail(dimv_);
}


inline void BackwardCorrection::backwardCorrectionParallel(
    const Robot& robot, SplitDirection& d, SplitSolution& s_new) const {
  d.split_direction.tail(dimKKT_-dimx_).noalias()
      = kkt_matrix_inverse_.bottomRightCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.u.noalias() -= d.du();
  robot.integrateConfiguration(d.dq(), -1, s_new.q);
  s_new.v.noalias() -= d.dv();
}


template <typename SplitSolutionType>
inline void BackwardCorrection::forwardCorrectionSerial(
    const Robot& robot, const SplitSolutionType& s_prev, 
    const SplitSolutionType& s_new_prev, SplitSolution& s_new) {
  robot.subtractConfiguration(s_new_prev.q, s_prev.q, 
                              x_res_.head(dimv_));
  x_res_.tail(dimv_) = s_new_prev.v - s_prev.v;
  dx_.noalias() = kkt_matrix_inverse_.bottomLeftCorner(dimx_, dimx_) * x_res_;
  robot.integrateConfiguration(dx_.head(dimv_), -1, s_new.q);
  s_new.v.noalias() -= dx_.tail(dimv_);
}


inline void BackwardCorrection::forwardCorrectionParallel(
    SplitDirection& d, SplitSolution& s_new) const {
  d.split_direction.head(dimKKT_-dimx_).noalias()
      = kkt_matrix_inverse_.topLeftCorner(dimKKT_-dimx_, dimx_) * x_res_;
  s_new.lmd.noalias() -= d.dlmd();
  s_new.gmm.noalias() -= d.dgmm();
  s_new.u.noalias() -= d.du();
}


inline void BackwardCorrection::computeDirection(const Robot& robot,
                                                 const SplitSolution& s, 
                                                 const SplitSolution& s_new, 
                                                 SplitDirection& d) {
  d.dlmd() = s_new.lmd - s.lmd;
  d.dgmm() = s_new.gmm - s.gmm;
  d.du() = s_new.u - s.u;
  robot.subtractConfiguration(s_new.q, s.q, d.dq());
  d.dv() = s_new.v - s.v;
}


template <typename MatrixType>
inline void BackwardCorrection::getAuxiliaryMatrix(
    const Eigen::MatrixBase<MatrixType>& aux_mat) const {
  assert(aux_mat.rows() == dimx_);
  assert(aux_mat.cols() == dimx_);
  (const_cast<Eigen::MatrixBase<MatrixType>&>(aux_mat)) 
      = - kkt_matrix_inverse_.topLeftCorner(dimx_, dimx_);
}

} // namespace idocp

#endif // IDOCP_BACKWARD_CORRECTION_HXX_ 