#include "idocp/cost/time_varying_configuration_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

TimeVaryingConfigurationCost::TimeVaryingConfigurationCost(const Robot& robot) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    t0_(0),
    q0_(Eigen::VectorXd::Zero(robot.dimq())),
    v0_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
}


TimeVaryingConfigurationCost::TimeVaryingConfigurationCost()
  : dimq_(0),
    dimv_(0),
    t0_(0),
    q0_(),
    v0_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    qf_weight_(),
    vf_weight_() {
}


TimeVaryingConfigurationCost::~TimeVaryingConfigurationCost() {
}


bool TimeVaryingConfigurationCost::useKinematics() const {
  return false;
}


void TimeVaryingConfigurationCost::set_ref(const double t0, 
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


void TimeVaryingConfigurationCost::set_q_weight(
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


void TimeVaryingConfigurationCost::set_v_weight(
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


void TimeVaryingConfigurationCost::set_a_weight(
    const Eigen::VectorXd& a_weight) {
  try {
    if (a_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: a_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  a_weight_ = a_weight;
}


void TimeVaryingConfigurationCost::set_qf_weight(
    const Eigen::VectorXd& qf_weight) {
  try {
    if (qf_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: qf_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  qf_weight_ = qf_weight;
}


void TimeVaryingConfigurationCost::set_vf_weight(
    const Eigen::VectorXd& vf_weight) {
  try {
    if (vf_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: vf_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  vf_weight_ = vf_weight;
}


double TimeVaryingConfigurationCost::l(Robot& robot, CostFunctionData& data, 
                                       const double t, const double dtau, 
                                       const SplitSolution& s) const {
  double l = 0;
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (q_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    l += (q_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  l += (v_weight_.array()*(s.v-v0_).array()*(s.v-v0_).array()).sum();
  l += (a_weight_.array()*s.a.array()*s.a.array()).sum();
  return 0.5 * dtau * l;
}


double TimeVaryingConfigurationCost::phi(Robot& robot, CostFunctionData& data, 
                                         const double t, 
                                         const SplitSolution& s) const {
  double phi = 0;
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    phi += (qf_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    phi += (qf_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  phi += (vf_weight_.array()*(s.v-v0_).array()*(s.v-v0_).array()).sum();
  return 0.5 * phi;
}


void TimeVaryingConfigurationCost::lq(Robot& robot, CostFunctionData& data, 
                                      const double t, const double dtau, 
                                      const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual) const {
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += dtau * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += dtau * q_weight_.array() * (s.q.array()-data.q_ref.array());
  }
}


void TimeVaryingConfigurationCost::lv(Robot& robot, CostFunctionData& data, 
                                      const double t, const double dtau, 
                                      const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += dtau * v_weight_.array() * (s.v.array()-v0_.array());
}


void TimeVaryingConfigurationCost::la(Robot& robot, CostFunctionData& data, 
                                      const double t, const double dtau, 
                                      const SplitSolution& s, 
                                      SplitKKTResidual& kkt_residual) const {
  kkt_residual.la.array() += dtau * a_weight_.array() * s.a.array();
}


void TimeVaryingConfigurationCost::lqq(Robot& robot, CostFunctionData& data, 
                                       const double t, const double dtau, 
                                       const SplitSolution& s, 
                                       SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += dtau * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += dtau * q_weight_;
  }
}


void TimeVaryingConfigurationCost::lvv(Robot& robot, CostFunctionData& data, 
                                       const double t, const double dtau, 
                                       const SplitSolution& s, 
                                       SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += dtau * v_weight_;
}


void TimeVaryingConfigurationCost::laa(Robot& robot, CostFunctionData& data, 
                                       const double t, const double dtau, 
                                       const SplitSolution& s, 
                                       SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qaa().diagonal().noalias() += dtau * a_weight_;
}


void TimeVaryingConfigurationCost::phiq(Robot& robot, CostFunctionData& data, 
                                        const double t, const SplitSolution& s, 
                                        SplitKKTResidual& kkt_residual) const {
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += qf_weight_.array() * (s.q.array()-data.q_ref.array());
  }
}


void TimeVaryingConfigurationCost::phiv(Robot& robot, CostFunctionData& data, 
                                        const double t, const SplitSolution& s, 
                                        SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += vf_weight_.array() * (s.v.array()-v0_.array());
}


void TimeVaryingConfigurationCost::phiqq(Robot& robot, CostFunctionData& data, 
                                         const double t, const SplitSolution& s, 
                                         SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
  }
}


void TimeVaryingConfigurationCost::phivv(Robot& robot, CostFunctionData& data, 
                                         const double t, const SplitSolution& s, 
                                         SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += vf_weight_;
}

} // namespace idocp