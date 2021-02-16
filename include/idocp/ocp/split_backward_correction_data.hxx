#ifndef IDOCP_SPLIT_BACKWARD_CORRECTION_DATA_HXX_
#define IDOCP_SPLIT_BACKWARD_CORRECTION_DATA_HXX_

#include "idocp/impulse/impulse_split_backward_correction_data.hpp"

#include <cassert>

namespace idocp {

inline SplitBackwardCorrectionData::SplitBackwardCorrectionData(
      const Robot& robot) 
  : dimv_(robot.dimv()),
    dimx_(2*robot.dimv()),
    dimu_(robot.dimu()),
    dimi_(0),
    dimKKT_(4*robot.dimv()+robot.dimu()),
    du_begin_(2*robot.dimv()), 
    dq_begin_(2*robot.dimv()+robot.dimu()), 
    dv_begin_(3*robot.dimv()+robot.dimu()),
    KKT_mat_inv_full_(Eigen::MatrixXd::Zero(
        4*robot.dimv()+robot.dimu()+robot.max_dimf(), 
        4*robot.dimv()+robot.dimu()+robot.max_dimf())),
    split_direction_full_(Eigen::VectorXd::Zero(
        4*robot.dimv()+robot.dimu()+robot.max_dimf())) {
}


inline SplitBackwardCorrectionData::SplitBackwardCorrectionData() 
  : dimv_(0),
    dimx_(0),
    dimu_(0),
    dimi_(0),
    dimKKT_(0),
    du_begin_(0), 
    dq_begin_(0), 
    dv_begin_(0),
    KKT_mat_inv_full_(),
    split_direction_full_() {
}


inline SplitBackwardCorrectionData::~SplitBackwardCorrectionData() {
}


inline void SplitBackwardCorrectionData::setImpulseStatus(
    const int dimi) {
  dimi_   = dimi;
  dimKKT_ = 2*dimx_ + dimu_ + dimi_;
  du_begin_ = dimx_ + dimi_;
  dq_begin_ = dimx_ + dimi_ + dimu_;
  dv_begin_ = dimx_ + dimi_ + dimu_ + dimv_;
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitBackwardCorrectionData::KKT_mat_inv() {
  return KKT_mat_inv_full_.topLeftCorner(dimKKT_, dimKKT_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitBackwardCorrectionData::KKT_mat_inv() const {
  return KKT_mat_inv_full_.topLeftCorner(dimKKT_, dimKKT_);
}


inline Eigen::Block<Eigen::MatrixXd> 
SplitBackwardCorrectionData::auxMat() {
  return KKT_mat_inv_full_.topLeftCorner(dimx_, dimx_);
}


inline const Eigen::Block<const Eigen::MatrixXd> 
SplitBackwardCorrectionData::auxMat() const {
  return KKT_mat_inv_full_.topLeftCorner(dimx_, dimx_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitBackwardCorrectionData::splitDirection() {
  return split_direction_full_.head(dimKKT_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitBackwardCorrectionData::splitDirection() const {
  return split_direction_full_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitBackwardCorrectionData::dlmd() {
  return split_direction_full_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitBackwardCorrectionData::dlmd() const {
  return split_direction_full_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitBackwardCorrectionData::dgmm() {
  return split_direction_full_.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitBackwardCorrectionData::dgmm() const {
  return split_direction_full_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitBackwardCorrectionData::dxi() {
  return split_direction_full_.segment(dimx_, dimi_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitBackwardCorrectionData::dxi() const {
  return split_direction_full_.segment(dimx_, dimi_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitBackwardCorrectionData::du() {
  return split_direction_full_.segment(du_begin_, dimu_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitBackwardCorrectionData::du() const {
  return split_direction_full_.segment(du_begin_, dimu_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitBackwardCorrectionData::dq() {
  return split_direction_full_.segment(dq_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitBackwardCorrectionData::dq() const {
  return split_direction_full_.segment(dq_begin_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> 
SplitBackwardCorrectionData::dv() {
  return split_direction_full_.segment(dv_begin_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitBackwardCorrectionData::dv() const {
  return split_direction_full_.segment(dv_begin_, dimv_);
}

} // namespace idocp

#endif // IDOCP_SPLIT_BACKWARD_CORRECTION_DATA_HXX_ 