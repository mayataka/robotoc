#ifndef IDOCP_SPLIT_DIRECTION_HPP_
#define IDOCP_SPLIT_DIRECTION_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class SplitDirection {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  SplitDirection(const Robot& robot) 
    : kkt_composition_(robot),
      split_direction_(Eigen::VectorXd::Zero(kkt_composition_.max_dimKKT())) {
  }

  SplitDirection() 
    : kkt_composition_(),
      split_direction_() {
  }

  ~SplitDirection() {
  }

  SplitDirection(const SplitDirection&) = default;

  SplitDirection& operator=(const SplitDirection&) = default;
 
  SplitDirection(SplitDirection&&) noexcept = default;

  SplitDirection& operator=(SplitDirection&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    kkt_composition_.set(robot);
  }

  inline Eigen::Ref<Eigen::VectorXd> split_direction() {
    return split_direction_.head(kkt_composition_.dimKKT());
  }

  inline Eigen::Ref<Eigen::VectorXd> dlmd() {
    return split_direction_.segment(kkt_composition_.Fq_begin(), 
                                    kkt_composition_.Fq_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dgmm() {
    return split_direction_.segment(kkt_composition_.Fv_begin(), 
                                    kkt_composition_.Fv_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dmu() {
    return split_direction_.segment(kkt_composition_.C_begin(), 
                                    kkt_composition_.C_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> da() {
    return split_direction_.segment(kkt_composition_.Qa_begin(), 
                                    kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> df() {
    return split_direction_.segment(kkt_composition_.Qf_begin(), 
                                    kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dq() {
    return split_direction_.segment(kkt_composition_.Qq_begin(), 
                                    kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dv() {
    return split_direction_.segment(kkt_composition_.Qv_begin(), 
                                    kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dx() {
    return split_direction_.segment(kkt_composition_.Qx_begin(), 
                                    kkt_composition_.Qx_size());
  }

  inline void setZero() {
    split_direction_.setZero();
  }

  inline int dimKKT() const {
    return kkt_composition_.dimKKT();
  }

  inline int max_dimKKT() const {
    return kkt_composition_.max_dimKKT();
  }

private:
  KKTComposition kkt_composition_;
  Eigen::VectorXd split_direction_;

};

} // namespace idocp 


#endif // IDOCP_SPLIT_OCP_DIRECTION_HPP_