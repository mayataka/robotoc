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
      split_direction_(Eigen::VectorXd::Zero(kkt_composition_.max_dimKKT())),
      du(robot.dimv()),
      dbeta(robot.dimv()),
      dimc_(robot.dim_passive()+robot.dimf()),
      dimf_(robot.dimf()) {
  }

  SplitDirection() 
    : kkt_composition_(),
      split_direction_(),
      du(),
      dbeta(),
      dimc_(0),
      dimf_(0) {
  }

  ~SplitDirection() {
  }

  SplitDirection(const SplitDirection&) = default;

  SplitDirection& operator=(const SplitDirection&) = default;
 
  SplitDirection(SplitDirection&&) noexcept = default;

  SplitDirection& operator=(SplitDirection&&) noexcept = default;

  inline void set(const Robot& robot) {
    kkt_composition_.setContactStatus(robot);
    dimc_ = robot.dim_passive() + robot.dimf();
    dimf_ = robot.dimf();
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

  inline Eigen::Ref<const Eigen::VectorXd> split_direction() const {
    return split_direction_.head(kkt_composition_.dimKKT());
  }

  inline Eigen::Ref<const Eigen::VectorXd> dlmd() const {
    return split_direction_.segment(kkt_composition_.Fq_begin(), 
                                    kkt_composition_.Fq_size());
  }

  inline Eigen::Ref<const Eigen::VectorXd> dgmm() const {
    return split_direction_.segment(kkt_composition_.Fv_begin(), 
                                    kkt_composition_.Fv_size());
  }

  inline Eigen::Ref<const Eigen::VectorXd> dmu() const {
    return split_direction_.segment(kkt_composition_.C_begin(), 
                                    kkt_composition_.C_size());
  }

  inline Eigen::Ref<const Eigen::VectorXd> da() const {
    return split_direction_.segment(kkt_composition_.Qa_begin(), 
                                    kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<const Eigen::VectorXd> df() const {
    return split_direction_.segment(kkt_composition_.Qf_begin(), 
                                    kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<const Eigen::VectorXd> dq() const {
    return split_direction_.segment(kkt_composition_.Qq_begin(), 
                                    kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<const Eigen::VectorXd> dv() const {
    return split_direction_.segment(kkt_composition_.Qv_begin(), 
                                    kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<const Eigen::VectorXd> dx() const {
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

  int dimc() const {
    return dimc_;
  }

  int dimf() const {
    return dimf_;
  }

  // condensed directions
  Eigen::VectorXd du, dbeta;

private:
  KKTComposition kkt_composition_;
  Eigen::VectorXd split_direction_;
  int dimc_, dimf_;

};

} // namespace idocp 


#endif // IDOCP_SPLIT_OCP_DIRECTION_HPP_