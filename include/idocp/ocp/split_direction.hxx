#ifndef IDOCP_SPLIT_DIRECTION_HXX_
#define IDOCP_SPLIT_DIRECTION_HXX_

namespace idocp {

inline SplitDirection::SplitDirection(const Robot& robot) 
  : du(robot.dimv()),
    dbeta(robot.dimv()),
    split_direction_(Eigen::VectorXd::Zero(
        5*robot.dimv()+robot.dim_passive()+2*robot.max_dimf())),
    dimv_(robot.dimv()), 
    dimx_(2*robot.dimv()), 
    dim_passive_(robot.dim_passive()), 
    max_dimf_(robot.max_dimf()), 
    dimf_(robot.dimf()), 
    max_dimc_(robot.dim_passive()+robot.max_dimf()), 
    dimc_(robot.dim_passive()+robot.dimf()),
    max_dimKKT_(5*robot.dimv()+robot.dim_passive()+2*robot.max_dimf()),
    dimKKT_(5*robot.dimv()+robot.dim_passive()+2*robot.dimf()) {
}


inline SplitDirection::SplitDirection() 
  : du(),
    dbeta(),
    split_direction_(),
    dimv_(0), 
    dimx_(0), 
    dim_passive_(0), 
    max_dimf_(0), 
    dimf_(0), 
    max_dimc_(0), 
    dimc_(0),
    max_dimKKT_(0),
    dimKKT_(0) {
}


inline SplitDirection::~SplitDirection() {
}


inline void SplitDirection::setContactStatus(const Robot& robot) {
  dimf_ = robot.dimf();
  dimc_ = dim_passive_ + robot.dimf();
  dimKKT_ = 5*dimv_ + dim_passive_ + 2*robot.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::split_direction() {
  return split_direction_.head(dimKKT_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dlmd() {
  return split_direction_.head(dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dgmm() {
  return split_direction_.segment(dimv_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dmu() {
  return split_direction_.segment(dimx_, dimc_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::da() {
  return split_direction_.segment(dimx_+dimc_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::df() {
  return split_direction_.segment(dimx_+dimc_+dimv_, dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dq() {
  return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dv() {
  return split_direction_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitDirection::dx() {
  return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::split_direction() const {
  return split_direction_.head(dimKKT_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dlmd() const {
  return split_direction_.head(dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dgmm() const {
  return split_direction_.segment(dimv_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dmu() const {
  return split_direction_.segment(dimx_, dimc_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::da() const {
  return split_direction_.segment(dimx_+dimc_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::df() const {
  return split_direction_.segment(dimx_+dimc_+dimv_, dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dq() const {
  return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dv() const {
  return split_direction_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitDirection::dx() const {
  return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
}


inline void SplitDirection::setZero() {
  split_direction_.setZero();
}


inline int SplitDirection::dimKKT() const {
  return dimKKT_;
}


inline int SplitDirection::max_dimKKT() const {
  return max_dimKKT_;
}


inline int SplitDirection::dimc() const {
  return dimc_;
}


inline int SplitDirection::dimf() const {
  return dimf_;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_OCP_DIRECTION_HXX_