#include "idocp/cost/time_varying_configuration_space_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

TimeVaryingConfigurationSpaceCost::TimeVaryingConfigurationSpaceCost(
    const Robot& robot) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    t0_(0),
    q0_(Eigen::VectorXd::Zero(robot.dimq())),
    v0_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qi_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vi_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    dvi_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
}


TimeVaryingConfigurationSpaceCost::TimeVaryingConfigurationSpaceCost()
  : dimq_(0),
    dimv_(0),
    t0_(0),
    q0_(),
    v0_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    qf_weight_(),
    vf_weight_(),
    qi_weight_(),
    vi_weight_(),
    dvi_weight_() {
}


TimeVaryingConfigurationSpaceCost::~TimeVaryingConfigurationSpaceCost() {
}


bool TimeVaryingConfigurationSpaceCost::useKinematics() const {
  return false;
}


void TimeVaryingConfigurationSpaceCost::set_ref(const double t0, 
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


void TimeVaryingConfigurationSpaceCost::set_q_weight(
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


void TimeVaryingConfigurationSpaceCost::set_v_weight(
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


void TimeVaryingConfigurationSpaceCost::set_a_weight(
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


void TimeVaryingConfigurationSpaceCost::set_qf_weight(
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


void TimeVaryingConfigurationSpaceCost::set_vf_weight(
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


void TimeVaryingConfigurationSpaceCost::set_qi_weight(
    const Eigen::VectorXd& qi_weight) {
  try {
    if (qi_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: qi_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  qi_weight_ = qi_weight;
}


void TimeVaryingConfigurationSpaceCost::set_vi_weight(
    const Eigen::VectorXd& vi_weight) {
  try {
    if (vi_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: vi_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  vi_weight_ = vi_weight;
}


void TimeVaryingConfigurationSpaceCost::set_dvi_weight(
    const Eigen::VectorXd& dvi_weight) {
  try {
    if (dvi_weight.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: dvi_weight.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  dvi_weight_ = dvi_weight;
}


double TimeVaryingConfigurationSpaceCost::computeStageCost(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
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


double TimeVaryingConfigurationSpaceCost::computeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
  double l = 0;
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (qf_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    l += (qf_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  l += (vf_weight_.array()*(s.v-v0_).array()*(s.v-v0_).array()).sum();
  return 0.5 * l;
}


double TimeVaryingConfigurationSpaceCost::computeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (qi_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    l += (qi_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  l += (vi_weight_.array()*(s.v-v0_).array()*(s.v-v0_).array()).sum();
  l += (dvi_weight_.array()*s.dv.array()*s.dv.array()).sum();
  return 0.5 * l;
}


void TimeVaryingConfigurationSpaceCost::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
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
  kkt_residual.lv().array()
      += dtau * v_weight_.array() * (s.v.array()-v0_.array());
  kkt_residual.la.array() += dtau * a_weight_.array() * s.a.array();
}


void TimeVaryingConfigurationSpaceCost::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
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
  kkt_residual.lv().array()
      += vf_weight_.array() * (s.v.array()-v0_.array());
}


void TimeVaryingConfigurationSpaceCost::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * qi_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += qi_weight_.array() * (s.q.array()-data.q_ref.array());
  }
  kkt_residual.lv().array()
      += vi_weight_.array() * (s.v.array()-v0_.array());
  kkt_residual.ldv.array() += dvi_weight_.array() * s.dv.array();
}


void TimeVaryingConfigurationSpaceCost::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += dtau * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += dtau * q_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += dtau * v_weight_;
  kkt_matrix.Qaa().diagonal().noalias() += dtau * a_weight_;
}


void TimeVaryingConfigurationSpaceCost::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += vf_weight_;
}


void TimeVaryingConfigurationSpaceCost::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.integrateConfiguration(q0_, v0_, t-t0_, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qi_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qi_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += vi_weight_;
  kkt_matrix.Qdvdv().diagonal().noalias() += dvi_weight_;
}

} // namespace idocp