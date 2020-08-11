#ifndef IDOCP_KKT_DIRECTION_HPP_
#define IDOCP_KKT_DIRECTION_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class KKTDirection {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  KKTDirection(const Robot& robot) 
    : kkt_composition_(robot),
      kkt_direction_(Eigen::VectorXd::Zero(kkt_composition_.max_dimKKT())) {
  }

  KKTDirection() 
    : kkt_composition_(),
      kkt_direction_() {
  }

  ~KKTDirection() {
  }

  KKTDirection(const KKTDirection&) = default;

  KKTDirection& operator=(const KKTDirection&) = default;
 
  KKTDirection(KKTDirection&&) noexcept = default;

  KKTDirection& operator=(KKTDirection&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    kkt_composition_.set(robot);
  }

  inline Eigen::Ref<Eigen::VectorXd> KKT_direction() {
    return kkt_direction_.head(kkt_composition_.dimKKT());
  }

  inline Eigen::Ref<Eigen::VectorXd> dlmd() {
    return kkt_direction_.segment(kkt_composition_.Fq_begin(), 
                                  kkt_composition_.Fq_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dgmm() {
    return kkt_direction_.segment(kkt_composition_.Fv_begin(), 
                                  kkt_composition_.Fv_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dmu() {
    return kkt_direction_.segment(kkt_composition_.C_begin(), 
                                  kkt_composition_.C_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> da() {
    return kkt_direction_.segment(kkt_composition_.Qa_begin(), 
                                  kkt_composition_.Qa_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> df() {
    return kkt_direction_.segment(kkt_composition_.Qf_begin(), 
                                  kkt_composition_.Qf_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dq() {
    return kkt_direction_.segment(kkt_composition_.Qq_begin(), 
                                  kkt_composition_.Qq_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dv() {
    return kkt_direction_.segment(kkt_composition_.Qv_begin(), 
                                  kkt_composition_.Qv_size());
  }

  inline Eigen::Ref<Eigen::VectorXd> dx() {
    return kkt_direction_.segment(kkt_composition_.Qx_begin(), 
                                  kkt_composition_.Qx_size());
  }

  inline void setZero() {
    kkt_direction_.setZero();
  }

  inline int dimKKT() const {
    return kkt_composition_.dimKKT();
  }

  inline int max_dimKKT() const {
    return kkt_composition_.max_dimKKT();
  }

private:
  KKTComposition kkt_composition_;
  Eigen::VectorXd kkt_direction_;

};

} // namespace idocp 


#endif // IDOCP_KKT_DIRECTION_HPP_