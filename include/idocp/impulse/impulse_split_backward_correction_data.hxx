#ifndef IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_DATA_HXX_
#define IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_DATA_HXX_

#include "idocp/impulse/impulse_split_backward_correction_data.hpp"

#include <cassert>

namespace idocp {

inline ImpulseSplitBackwardCorrectionData::ImpulseSplitBackwardCorrectionData(
      const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimf_(0),
    dimKKT_(5*robot.dimv()),
    df_begin_(2*robot.dimv()), 
    dq_begin_(2*robot.dimv()), 
    dv_begin_(3*robot.dimv()),
    KKT_mat_inv_full_(Eigen::MatrixXd::Zero(4*robot.dimv()+2*robot.max_dimf(), 
                                            4*robot.dimv()+2*robot.max_dimf())),
    split_direction_full_(
        Eigen::VectorXd::Zero(4*robot.dimv()+2*robot.max_dimf())) {
}


inline ImpulseSplitBackwardCorrectionData::ImpulseSplitBackwardCorrectionData() 
  : dimv_(0),
    dimx_(0),
    dimf_(0),
    dimKKT_(0),
    df_begin_(0), 
    dq_begin_(0), 
    dv_begin_(0),
    KKT_mat_inv_full_(),
    split_direction_full_() {
}


inline ImpulseSplitBackwardCorrectionData::
~ImpulseSplitBackwardCorrectionData() {
}


inline void ImpulseSplitBackwardCorrectionData::setImpulseStatus(
    const int dimf) {
  dimf_   = dimf;
  dimKKT_ = 2*dimx_ + 2*dimf_;
  df_begin_ = dimx_ + dimf_;
  dq_begin_ = dimx_ + 2*dimf_;
  dv_begin_ = dimx_ + 2*dimf_ + dimv_;
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseSplitBackwardCorrectionData::KKT_mat_inv() {
  return KKT_mat_inv_full_.topLeftCorner(dimKKT_, dimKKT_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitBackwardCorrectionData::KKT_mat_inv() const {
  return KKT_mat_inv_full_.topLeftCorner(dimKKT_, dimKKT_);
}


inline Eigen::Block<Eigen::MatrixXd> 
ImpulseSplitBackwardCorrectionData::auxMat() {
  return KKT_mat_inv_full_.topLeftCorner(dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
ImpulseSplitBackwardCorrectionData::auxMat() const {
  return KKT_mat_inv_full_.topLeftCorner(dimx_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::splitDirection() {
  return split_direction_full_.head(dimKKT_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::splitDirection() const {
  return split_direction_full_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dlmd() {
  return split_direction_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dlmd() const {
  return split_direction_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dgmm() {
  return split_direction_full_.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dgmm() const {
  return split_direction_full_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dmu() {
  return split_direction_full_.segment(dimx_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dmu() const {
  return split_direction_full_.segment(dimx_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::df() {
  return split_direction_full_.segment(df_begin_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::df() const {
  return split_direction_full_.segment(df_begin_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dq() {
  return split_direction_full_.segment(dq_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dq() const {
  return split_direction_full_.segment(dq_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dv() {
  return split_direction_full_.segment(dv_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
ImpulseSplitBackwardCorrectionData::dv() const {
  return split_direction_full_.segment(dv_begin_, dimv_);
}

} // namespace idocp

#endif // IDOCP_IMPULSE_SPLIT_BACKWARD_CORRECTION_DATA_HXX_ 