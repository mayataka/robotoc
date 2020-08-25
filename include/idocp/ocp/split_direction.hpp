#ifndef IDOCP_SPLIT_DIRECTION_HPP_
#define IDOCP_SPLIT_DIRECTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class SplitDirection {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  SplitDirection(const Robot& robot) 
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

  SplitDirection() 
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

  ~SplitDirection() {
  }

  SplitDirection(const SplitDirection&) = default;

  SplitDirection& operator=(const SplitDirection&) = default;
 
  SplitDirection(SplitDirection&&) noexcept = default;

  SplitDirection& operator=(SplitDirection&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    dimf_ = robot.dimf();
    dimc_ = dim_passive_ + robot.dimf();
    dimKKT_ = 5*dimv_ + dim_passive_ + 2*robot.dimf();
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> split_direction() {
    return split_direction_.head(dimKKT_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> dlmd() {
    return split_direction_.head(dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> dgmm() {
    return split_direction_.segment(dimv_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> dmu() {
    return split_direction_.segment(dimx_, dimc_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> da() {
    return split_direction_.segment(dimx_+dimc_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> df() {
    return split_direction_.segment(dimx_+dimc_+dimv_, dimf_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> dq() {
    return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> dv() {
    return split_direction_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
  }

  inline Eigen::VectorBlock<Eigen::VectorXd> dx() {
    return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> 
  split_direction() const {
    return split_direction_.head(dimKKT_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> dlmd() const {
    return split_direction_.head(dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> dgmm() const {
    return split_direction_.segment(dimv_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> dmu() const {
    return split_direction_.segment(dimx_, dimc_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> da() const {
    return split_direction_.segment(dimx_+dimc_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> df() const {
    return split_direction_.segment(dimx_+dimc_+dimv_, dimf_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> dq() const {
    return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> dv() const {
    return split_direction_.segment(dimx_+dimc_+2*dimv_+dimf_, dimv_);
  }

  inline const Eigen::VectorBlock<const Eigen::VectorXd> dx() const {
    return split_direction_.segment(dimx_+dimc_+dimv_+dimf_, dimx_);
  }

  inline void setZero() {
    split_direction_.setZero();
  }

  inline int dimKKT() const {
    return dimKKT_;
  }

  inline int max_dimKKT() const {
    return max_dimKKT_;
  }

  inline int dimc() const {
    return dimc_;
  }

  inline int dimf() const {
    return dimf_;
  }

  // condensed directions
  Eigen::VectorXd du, dbeta;

private:
  Eigen::VectorXd split_direction_;
  int dimv_, dimx_, dim_passive_, max_dimf_, dimf_, max_dimc_, dimc_, dimKKT_, 
      max_dimKKT_;

};

} // namespace idocp 


#endif // IDOCP_SPLIT_OCP_DIRECTION_HPP_