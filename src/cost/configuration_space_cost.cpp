#include "idocp/cost/configuration_space_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

ConfigurationSpaceCost::ConfigurationSpaceCost(const Robot& robot)
  : CostFunctionComponentBase(),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    q_ref_(Eigen::VectorXd::Zero(robot.dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    u_ref_(Eigen::VectorXd::Zero(robot.dimu())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    u_weight_(Eigen::VectorXd::Zero(robot.dimu())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qi_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vi_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    dvi_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
  if (robot.hasFloatingBase()) {
    robot.normalizeConfiguration(q_ref_);
  }
}


ConfigurationSpaceCost::ConfigurationSpaceCost()
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


ConfigurationSpaceCost::~ConfigurationSpaceCost() {
}


bool ConfigurationSpaceCost::useKinematics() const {
  return false;
}


void ConfigurationSpaceCost::set_q_ref(const Eigen::VectorXd& q_ref) {
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
  q_ref_ = q_ref;
}


void ConfigurationSpaceCost::set_v_ref(const Eigen::VectorXd& v_ref) {
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
  v_ref_ = v_ref;
}


void ConfigurationSpaceCost::set_u_ref(const Eigen::VectorXd& u_ref) {
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
  u_ref_ = u_ref;
}


void ConfigurationSpaceCost::set_q_weight(const Eigen::VectorXd& q_weight) {
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


void ConfigurationSpaceCost::set_v_weight(const Eigen::VectorXd& v_weight) {
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


void ConfigurationSpaceCost::set_a_weight(const Eigen::VectorXd& a_weight) {
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


void ConfigurationSpaceCost::set_u_weight(const Eigen::VectorXd& u_weight) {
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
  u_weight_ = u_weight;
}


void ConfigurationSpaceCost::set_qf_weight(const Eigen::VectorXd& qf_weight) {
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


void ConfigurationSpaceCost::set_vf_weight(const Eigen::VectorXd& vf_weight) {
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


void ConfigurationSpaceCost::set_qi_weight(const Eigen::VectorXd& qi_weight) {
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


void ConfigurationSpaceCost::set_vi_weight(const Eigen::VectorXd& vi_weight) {
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


void ConfigurationSpaceCost::set_dvi_weight(const Eigen::VectorXd& dvi_weight) {
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


double ConfigurationSpaceCost::computeStageCost(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s) const {
  double l = 0;
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    l += (q_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  else {
    l += (q_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  l += (v_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  l += (a_weight_.array()*s.a.array()*s.a.array()).sum();
  l += (u_weight_.array()*(s.u-u_ref_).array()*(s.u-u_ref_).array()).sum();
  return 0.5 * dt * l;
}


double ConfigurationSpaceCost::computeTerminalCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s) const {
double l = 0;
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    l += (qf_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  else {
    l += (qf_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  l += (vf_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  return 0.5 * l;
}


double ConfigurationSpaceCost::computeImpulseCost(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s) const {
  double l = 0;
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    l += (qi_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  else {
    l += (qi_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  l += (vi_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  l += (dvi_weight_.array()*s.dv.array()*s.dv.array()).sum();
  return 0.5 * l;
}


void ConfigurationSpaceCost::computeStageCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_residual.lq().noalias()
        += dt * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += dt * q_weight_.array() * (s.q.array()-q_ref_.array());
  }
  kkt_residual.lv().array()
      += dt * v_weight_.array() * (s.v.array()-v_ref_.array());
  kkt_residual.la.array() += dt * a_weight_.array() * s.a.array();
  kkt_residual.lu().array() 
      += dt * u_weight_.array() * (s.u.array()-u_ref_.array());
}


void ConfigurationSpaceCost::computeTerminalCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual) const {
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += qf_weight_.array() * (s.q.array()-q_ref_.array());
  }
  kkt_residual.lv().array()
      += vf_weight_.array() * (s.v.array()-v_ref_.array());
}


void ConfigurationSpaceCost::computeImpulseCostDerivatives(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& 
    s, ImpulseSplitKKTResidual& kkt_residual) const {
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * qi_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += qi_weight_.array() * (s.q.array()-q_ref_.array());
  }
  kkt_residual.lv().array()
      += vi_weight_.array() * (s.v.array()-v_ref_.array());
  kkt_residual.ldv.array() += dvi_weight_.array() * s.dv.array();
}


void ConfigurationSpaceCost::computeStageCostHessian(
    Robot& robot, CostFunctionData& data, const double t, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += dt * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += dt * q_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += dt * v_weight_;
  kkt_matrix.Qaa().diagonal().noalias() += dt * a_weight_;
  kkt_matrix.Quu().diagonal().noalias() += dt * u_weight_;
}


void ConfigurationSpaceCost::computeTerminalCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
  }
  kkt_matrix.Qvv().diagonal().noalias() += vf_weight_;
}


void ConfigurationSpaceCost::computeImpulseCostHessian(
    Robot& robot, CostFunctionData& data, const double t, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
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