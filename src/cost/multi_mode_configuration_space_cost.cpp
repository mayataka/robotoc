#include "robotoc/cost/multi_mode_configuration_space_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace robotoc {

MultiModeConfigurationSpaceCost::MultiModeConfigurationSpaceCost(
    const Robot& robot)
  : CostFunctionComponentBase(),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    q_ref_(),
    v_ref_(),
    u_ref_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    u_weight_(),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qi_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vi_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    dvi_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
  q_ref_[0] = Eigen::VectorXd::Zero(robot.dimq());
  v_ref_[0] = Eigen::VectorXd::Zero(robot.dimv());
  u_ref_[0] = Eigen::VectorXd::Zero(robot.dimu());
  q_weight_[0] = Eigen::VectorXd::Zero(robot.dimv());
  v_weight_[0] = Eigen::VectorXd::Zero(robot.dimv());
  a_weight_[0] = Eigen::VectorXd::Zero(robot.dimv());
  u_weight_[0] = Eigen::VectorXd::Zero(robot.dimu());
  if (robot.hasFloatingBase()) {
    robot.normalizeConfiguration(q_ref_[0]);
  }
}


MultiModeConfigurationSpaceCost::MultiModeConfigurationSpaceCost()
  : CostFunctionComponentBase(),
    dimq_(0),
    dimv_(0),
    dimu_(),
    q_ref_(),
    v_ref_(),
    u_ref_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    u_weight_(),
    qf_weight_(),
    vf_weight_(),
    qi_weight_(),
    vi_weight_(),
    dvi_weight_() {
}


MultiModeConfigurationSpaceCost::~MultiModeConfigurationSpaceCost() {
}


void MultiModeConfigurationSpaceCost::set_q_ref(const Eigen::VectorXd& q_ref,
                                                const int contact_mode_id) {
  try {
    if (q_ref.size() != dimq_) {
      throw std::invalid_argument(
          "invalid size: q_ref.size() must be " + std::to_string(dimq_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  q_ref_[contact_mode_id] = q_ref;
}


void MultiModeConfigurationSpaceCost::set_q_ref(
    const Eigen::VectorXd& q_ref, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_q_ref(q_ref, e);
  }
}


void MultiModeConfigurationSpaceCost::set_v_ref(const Eigen::VectorXd& v_ref,
                                                const int contact_mode_id) {
  try {
    if (v_ref.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: v_ref.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  v_ref_[contact_mode_id] = v_ref;
}


void MultiModeConfigurationSpaceCost::set_v_ref(
    const Eigen::VectorXd& v_ref, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_v_ref(v_ref, e);
  }
}


void MultiModeConfigurationSpaceCost::set_u_ref(const Eigen::VectorXd& u_ref,
                                                const int contact_mode_id) {
  try {
    if (u_ref.size() != dimu_) {
      throw std::invalid_argument(
          "invalid size: u_ref.size() must be " + std::to_string(dimu_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  u_ref_[contact_mode_id] = u_ref;
}


void MultiModeConfigurationSpaceCost::set_u_ref(
    const Eigen::VectorXd& u_ref, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_u_ref(u_ref, e);
  }
}


void MultiModeConfigurationSpaceCost::set_q_weight(const Eigen::VectorXd& q_weight,
                                                   const int contact_mode_id) {
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
  q_weight_[contact_mode_id] = q_weight;
}


void MultiModeConfigurationSpaceCost::set_q_weight(
    const Eigen::VectorXd& q_weight, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_q_weight(q_weight, e);
  }
}


void MultiModeConfigurationSpaceCost::set_v_weight(const Eigen::VectorXd& v_weight,
                                                   const int contact_mode_id) {
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
  v_weight_[contact_mode_id] = v_weight;
}


void MultiModeConfigurationSpaceCost::set_v_weight(
    const Eigen::VectorXd& v_weight, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_v_weight(v_weight, e);
  }
}


void MultiModeConfigurationSpaceCost::set_a_weight(const Eigen::VectorXd& a_weight,
                                                   const int contact_mode_id) {
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
  a_weight_[contact_mode_id] = a_weight;
}


void MultiModeConfigurationSpaceCost::set_a_weight(
    const Eigen::VectorXd& a_weight, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_a_weight(a_weight, e);
  }
}


void MultiModeConfigurationSpaceCost::set_u_weight(const Eigen::VectorXd& u_weight,
                                                   const int contact_mode_id) {
  try {
    if (u_weight.size() != dimu_) {
      throw std::invalid_argument(
          "invalid size: u_weight.size() must be " + std::to_string(dimu_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  u_weight_[contact_mode_id] = u_weight;
}


void MultiModeConfigurationSpaceCost::set_u_weight(
    const Eigen::VectorXd& u_weight, const std::vector<int>& contact_mode_ids) {
  for (const int e : contact_mode_ids) {
    set_u_weight(u_weight, e);
  }
}


void MultiModeConfigurationSpaceCost::set_q_weight_terminal(const Eigen::VectorXd& q_weight_terminal) {
  try {
    if (q_weight_terminal.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: q_weight_terminal.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  qf_weight_ = q_weight_terminal;
}


void MultiModeConfigurationSpaceCost::set_v_weight_terminal(const Eigen::VectorXd& v_weight_terminal) {
  try {
    if (v_weight_terminal.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: v_weight_terminal.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  vf_weight_ = v_weight_terminal;
}


void MultiModeConfigurationSpaceCost::set_q_weight_impulse(const Eigen::VectorXd& q_weight_impulse) {
  try {
    if (q_weight_impulse.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: q_weight_impulse.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  qi_weight_ = q_weight_impulse;
}


void MultiModeConfigurationSpaceCost::set_v_weight_impulse(const Eigen::VectorXd& v_weight_impulse) {
  try {
    if (v_weight_impulse.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: v_weight_impulse.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  vi_weight_ = v_weight_impulse;
}


void MultiModeConfigurationSpaceCost::set_dv_weight_impulse(const Eigen::VectorXd& dv_weight_impulse) {
  try {
    if (dv_weight_impulse.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: dv_weight_impulse.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  dvi_weight_ = dv_weight_impulse;
}


bool MultiModeConfigurationSpaceCost::useKinematics() const {
  return false;
}


double MultiModeConfigurationSpaceCost::evalStageCost(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& q_ref = q_ref_.at(contact_mode_id);
  const auto& v_ref = v_ref_.at(contact_mode_id);
  const auto& u_ref = u_ref_.at(contact_mode_id);
  const auto& q_weight = q_weight_.at(contact_mode_id);
  const auto& v_weight = v_weight_.at(contact_mode_id);
  const auto& a_weight = a_weight_.at(contact_mode_id);
  const auto& u_weight = u_weight_.at(contact_mode_id);
  double l = 0;
  robot.subtractConfiguration(s.q, q_ref, data.qdiff);
  l += (q_weight.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  l += (v_weight.array()*(s.v-v_ref).array()*(s.v-v_ref).array()).sum();
  l += (a_weight.array()*s.a.array()*s.a.array()).sum();
  l += (u_weight.array()*(s.u-u_ref).array()*(s.u-u_ref).array()).sum();
  return 0.5 * grid_info.dt * l;
}


void MultiModeConfigurationSpaceCost::evalStageCostDerivatives(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTResidual& kkt_residual) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& q_ref = q_ref_.at(contact_mode_id);
  const auto& v_ref = v_ref_.at(contact_mode_id);
  const auto& u_ref = u_ref_.at(contact_mode_id);
  const auto& q_weight = q_weight_.at(contact_mode_id);
  const auto& v_weight = v_weight_.at(contact_mode_id);
  const auto& a_weight = a_weight_.at(contact_mode_id);
  const auto& u_weight = u_weight_.at(contact_mode_id);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += grid_info.dt * data.J_qdiff.transpose() * q_weight.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array() += grid_info.dt * q_weight.array() * data.qdiff.array();
  }
  kkt_residual.lv().array()
      += grid_info.dt * v_weight.array() * (s.v.array()-v_ref.array());
  kkt_residual.la.array() += grid_info.dt * a_weight.array() * s.a.array();
  kkt_residual.lu.array() 
      += grid_info.dt * u_weight.array() * (s.u.array()-u_ref.array());
}


void MultiModeConfigurationSpaceCost::evalStageCostHessian(
    Robot& robot, const ContactStatus& contact_status, CostFunctionData& data, 
    const GridInfo& grid_info, const SplitSolution& s, 
    SplitKKTMatrix& kkt_matrix) const {
  const int contact_mode_id = contact_status.contactModeId();
  const auto& q_weight = q_weight_.at(contact_mode_id);
  const auto& v_weight = v_weight_.at(contact_mode_id);
  const auto& a_weight = a_weight_.at(contact_mode_id);
  const auto& u_weight = u_weight_.at(contact_mode_id);
  if (robot.hasFloatingBase()) {
    kkt_matrix.Qqq().noalias()
        += grid_info.dt * data.J_qdiff.transpose() * q_weight.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += grid_info.dt * q_weight;
  }
  kkt_matrix.Qvv().diagonal().noalias() += grid_info.dt * v_weight;
  kkt_matrix.Qaa.diagonal().noalias() += grid_info.dt * a_weight;
  kkt_matrix.Quu.diagonal().noalias() += grid_info.dt * u_weight;
}


double MultiModeConfigurationSpaceCost::evalTerminalCost(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s) const {
  const auto& q_ref = q_ref_.at(0);
  const auto& v_ref = v_ref_.at(0);
  double l = 0;
  robot.subtractConfiguration(s.q, q_ref, data.qdiff);
  l += (qf_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  l += (vf_weight_.array()*(s.v-v_ref).array()*(s.v-v_ref).array()).sum();
  return 0.5 * l;
}


void MultiModeConfigurationSpaceCost::evalTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  const auto& q_ref = q_ref_.at(0);
  const auto& v_ref = v_ref_.at(0);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array() += qf_weight_.array() * data.qdiff.array();
  }
  kkt_residual.lv().array()
      += vf_weight_.array() * (s.v.array()-v_ref.array());
}


void MultiModeConfigurationSpaceCost::evalTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const GridInfo& grid_info, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += vf_weight_;
}


double MultiModeConfigurationSpaceCost::evalImpulseCost(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s) const {
  const int impulse_mode_id = impulse_status.impulseModeId();
  const auto& q_ref = q_ref_.at(impulse_mode_id);
  const auto& v_ref = v_ref_.at(impulse_mode_id);
  double l = 0;
  robot.subtractConfiguration(s.q, q_ref, data.qdiff);
  l += (qi_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  l += (vi_weight_.array()*(s.v-v_ref).array()*(s.v-v_ref).array()).sum();
  l += (dvi_weight_.array()*s.dv.array()*s.dv.array()).sum();
  return 0.5 * l;
}


void MultiModeConfigurationSpaceCost::evalImpulseCostDerivatives(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTResidual& kkt_residual) const {
  const int impulse_mode_id = impulse_status.impulseModeId();
  const auto& q_ref = q_ref_.at(impulse_mode_id);
  const auto& v_ref = v_ref_.at(impulse_mode_id);
  if (robot.hasFloatingBase()) {
    robot.dSubtractConfiguration_dqf(s.q, q_ref, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * qi_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array() += qi_weight_.array() * data.qdiff.array();
  }
  kkt_residual.lv().array()
      += vi_weight_.array() * (s.v.array()-v_ref.array());
  kkt_residual.ldv.array() += dvi_weight_.array() * s.dv.array();
}


void MultiModeConfigurationSpaceCost::evalImpulseCostHessian(
    Robot& robot, const ImpulseStatus& impulse_status, CostFunctionData& data, 
    const GridInfo& grid_info, const ImpulseSplitSolution& s, 
    ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qi_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qi_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += vi_weight_;
  kkt_matrix.Qdvdv.diagonal().noalias() += dvi_weight_;
}

} // namespace robotoc