#include "idocp/cost/trotting_configuration_space_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

TrottingConfigurationSpaceCost::TrottingConfigurationSpaceCost(
    const Robot& robot) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    t_start_(0), 
    t_period_(0), 
    step_length_(0), 
    v_com_(0),
    q_standing_(Eigen::VectorXd::Zero(robot.dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())), 
    swing_angles_() {
}


TrottingConfigurationSpaceCost::TrottingConfigurationSpaceCost()
  : dimq_(0),
    dimv_(0),
    t_start_(0), 
    t_period_(0), 
    step_length_(0), 
    v_com_(0),
    q_standing_(),
    v_ref_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    qf_weight_(),
    vf_weight_(),
    swing_angles_() {
}


TrottingConfigurationSpaceCost::~TrottingConfigurationSpaceCost() {
}


bool TrottingConfigurationSpaceCost::useKinematics() const {
  return false;
}


void TrottingConfigurationSpaceCost::set_ref(
    const double t_start, const double t_period, 
    const Eigen::VectorXd& q_standing, const double step_length,
    const TrottingSwingAngles& swing_angles) {
  try {
    if (q_standing.size() != dimq_) {
      throw std::invalid_argument(
          "invalid size: q_standing.size() must be " + std::to_string(dimq_) + "!");
    }
    if (t_period <= 0) {
      throw std::invalid_argument(
          "invalid argment: t_period must be positive!");
    }
    if (step_length <= 0) {
      throw std::invalid_argument(
          "invalid argment: step_length must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  t_start_ = t_start;
  t_period_ = t_period;
  q_standing_ = q_standing;
  step_length_ = step_length;
  v_com_ = step_length / t_period;
  v_ref_.setZero();
  v_ref_.coeffRef(0) = v_com_;
  swing_angles_ = swing_angles;
}


void TrottingConfigurationSpaceCost::set_q_weight(
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


void TrottingConfigurationSpaceCost::set_v_weight(
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


void TrottingConfigurationSpaceCost::set_a_weight(
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


void TrottingConfigurationSpaceCost::set_qf_weight(
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


void TrottingConfigurationSpaceCost::set_vf_weight(
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


void TrottingConfigurationSpaceCost::set_qi_weight(
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


void TrottingConfigurationSpaceCost::set_vi_weight(
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


void TrottingConfigurationSpaceCost::set_dvi_weight(
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


double TrottingConfigurationSpaceCost::computeStageCost(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s) const {
  double l = 0;
  update_q_ref(t, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (q_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    l += (q_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  l += (v_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  l += (a_weight_.array()*s.a.array()*s.a.array()).sum();
  return 0.5 * dt * l;
}


double TrottingConfigurationSpaceCost::computeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
  double l = 0;
  update_q_ref(t, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (qf_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    l += (qf_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  l += (vf_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  return 0.5 * l;
}


double TrottingConfigurationSpaceCost::computeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  update_q_ref(t, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    l += (qi_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    l += (qi_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  l += (vi_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  l += (dvi_weight_.array()*s.dv.array()*s.dv.array()).sum();
  return 0.5 * l;
}


void TrottingConfigurationSpaceCost::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  update_q_ref(t, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += dt * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += dt * q_weight_.array() * (s.q.array()-data.q_ref.array());
  }
  kkt_residual.lv().array()
      += dt * v_weight_.array() * (s.v.array()-v_ref_.array());
  kkt_residual.la.array() += dt * a_weight_.array() * s.a.array();
}


void TrottingConfigurationSpaceCost::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  update_q_ref(t, data.q_ref);
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
      += vf_weight_.array() * (s.v.array()-v_ref_.array());
}


void TrottingConfigurationSpaceCost::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  update_q_ref(t, data.q_ref);
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
      += vi_weight_.array() * (s.v.array()-v_ref_.array());
  kkt_residual.ldv.array() += dvi_weight_.array() * s.dv.array();
}


void TrottingConfigurationSpaceCost::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    update_q_ref(t, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += dt * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += dt * q_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += dt * v_weight_;
  kkt_matrix.Qaa().diagonal().noalias() += dt * a_weight_;
}


void TrottingConfigurationSpaceCost::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    update_q_ref(t, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += vf_weight_;
}


void TrottingConfigurationSpaceCost::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    update_q_ref(t, data.q_ref);
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