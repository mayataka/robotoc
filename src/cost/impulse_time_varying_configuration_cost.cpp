#include "idocp/cost/impulse_time_varying_configuration_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

ImpulseTimeVaryingConfigurationCost::ImpulseTimeVaryingConfigurationCost(
    const Robot& robot) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    t0_(0),
    q0_(Eigen::VectorXd::Zero(robot.dimq())),
    v0_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    dv_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
}


ImpulseTimeVaryingConfigurationCost::ImpulseTimeVaryingConfigurationCost()
  : dimq_(0),
    dimv_(0),
    t0_(0),
    q0_(),
    v0_(),
    q_weight_(),
    v_weight_(),
    dv_weight_() {
}


ImpulseTimeVaryingConfigurationCost::~ImpulseTimeVaryingConfigurationCost() {
}


void ImpulseTimeVaryingConfigurationCost::set_ref(const double t0, 
                                                  const Eigen::VectorXd q0, 
                                                  const Eigen::VectorXd v0) {
  t0_ = t0;
  if (q0.size() == dimq_) {
    q0_ = q0;
  }
  else {
    std::cout << "invalid argment in set_ref(): size of q0 must be " 
              << dimq_ << std::endl;
  }
  if (v0.size() == dimv_) {
    v0_ = v0;
  }
  else {
    std::cout << "invalid argment in set_ref(): size of v0 must be " 
              << dimv_ << std::endl;
  }
}


void ImpulseTimeVaryingConfigurationCost::set_q_weight(
    const Eigen::VectorXd& q_weight) {
  if (q_weight.size() == dimv_) {
    q_weight_ = q_weight;
  }
  else {
    std::cout << "invalid argment in set_q_weight(): size of q_weight must be " 
              << dimv_ << std::endl;
  }
}


void ImpulseTimeVaryingConfigurationCost::set_v_weight(
    const Eigen::VectorXd& v_weight) {
  if (v_weight.size() == dimv_) {
    v_weight_ = v_weight;
  }
  else {
    std::cout << "invalid argment in set_v_weight(): size of v_weight must be " 
              << dimv_ << std::endl;
  }
}


void ImpulseTimeVaryingConfigurationCost::set_dv_weight(
    const Eigen::VectorXd& dv_weight) {
  if (dv_weight.size() == dimv_) {
    dv_weight_ = dv_weight;
  }
  else {
    std::cout << "invalid argment in set_dv_weight(): size of dv_weight must be " 
              << dimv_ << std::endl;
  }
}


double ImpulseTimeVaryingConfigurationCost::l(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (q_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  else {
    l += (q_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  l += (v_weight_.array()*(s.v-v0_).array()*(s.v-v0_).array()).sum();
  l += (dv_weight_.array()*s.dv.array()*s.dv.array()).sum();
  return 0.5 * l;
}


void ImpulseTimeVaryingConfigurationCost::lq(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += q_weight_.array() * (s.q.array()-data.q_ref.array());
  }
}


void ImpulseTimeVaryingConfigurationCost::lv(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += v_weight_.array() * (s.v.array()-v0_.array());
}


void ImpulseTimeVaryingConfigurationCost::ldv(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  kkt_residual.ldv.array() += dv_weight_.array() * s.dv.array();
}


void ImpulseTimeVaryingConfigurationCost::lqq(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (robot.has_floating_base()) {
    robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += q_weight_;
  }
}


void ImpulseTimeVaryingConfigurationCost::lvv(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += v_weight_;
}


void ImpulseTimeVaryingConfigurationCost::ldvdv(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qdvdv().diagonal().noalias() += dv_weight_;
}

} // namespace idocp