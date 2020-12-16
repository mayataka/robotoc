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
  try {
    if (q0.size() != dimq_) {
      throw std::invalid_argument(
          "invalid size: q0.size() must be " + std::to_string(dimq_) + "!");
    }
    if (v0.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: v0.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  t0_ = t0;
  q0_ = q0;
  v0_ = v0;
}


void ImpulseTimeVaryingConfigurationCost::set_q_weight(
    const Eigen::VectorXd& q_weight) {
  try {
    if (q_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: q_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  q_weight_ = q_weight;
}


void ImpulseTimeVaryingConfigurationCost::set_v_weight(
    const Eigen::VectorXd& v_weight) {
  try {
    if (v_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: v_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  v_weight_ = v_weight;
}


void ImpulseTimeVaryingConfigurationCost::set_dv_weight(
    const Eigen::VectorXd& dv_weight) {
  try {
    if (dv_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: dv_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  dv_weight_ = dv_weight;
}


double ImpulseTimeVaryingConfigurationCost::l(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
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
  if (robot.hasFloatingBase()) {
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
  if (robot.hasFloatingBase()) {
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