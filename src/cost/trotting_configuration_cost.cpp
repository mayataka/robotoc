#include "idocp/cost/trotting_configuration_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

TrottingConfigurationCost::TrottingConfigurationCost(const Robot& robot) 
  : dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    t_start_(0), 
    t_period_(0), 
    step_length_(0), 
    v_com_(0),
    q_standing_(Eigen::VectorXd::Zero(robot.dimq())),
    q_even_step_(Eigen::VectorXd::Zero(robot.dimq())),
    q_odd_step_(Eigen::VectorXd::Zero(robot.dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
}


TrottingConfigurationCost::TrottingConfigurationCost()
  : dimq_(0),
    dimv_(0),
    t_start_(0), 
    t_period_(0), 
    step_length_(0), 
    v_com_(0),
    q_standing_(),
    q_even_step_(),
    q_odd_step_(),
    v_ref_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    qf_weight_(),
    vf_weight_() {
}


TrottingConfigurationCost::~TrottingConfigurationCost() {
}


bool TrottingConfigurationCost::useKinematics() const {
  return false;
}


void TrottingConfigurationCost::set_ref(const double t_start, 
                                        const double t_period, 
                                        const Eigen::VectorXd q_standing, 
                                        const double step_length,
                                        const double front_thigh_swing_angle, 
                                        const double front_knee_swing_angle, 
                                        const double hip_thigh_swing_angle, 
                                        const double hip_knee_swing_angle) {
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
  q_even_step_ = q_standing;
  q_odd_step_ = q_standing;

  // swing front leg
  q_even_step_.coeffRef(8) -= front_thigh_swing_angle;
  q_even_step_.coeffRef(9) -= front_knee_swing_angle;
  q_odd_step_.coeffRef(14) -= front_thigh_swing_angle;
  q_odd_step_.coeffRef(15) -= front_knee_swing_angle;

  // stance front leg
  q_even_step_.coeffRef(14) += front_thigh_swing_angle;
  q_even_step_.coeffRef(15) -= front_knee_swing_angle;
  q_odd_step_.coeffRef(8)   += front_thigh_swing_angle;
  q_odd_step_.coeffRef(9)   -= front_knee_swing_angle;

  // stance hip leg
  q_odd_step_.coeffRef(14) += hip_thigh_swing_angle;
  q_odd_step_.coeffRef(15) += hip_knee_swing_angle;
  q_odd_step_.coeffRef(17) += hip_thigh_swing_angle;
  q_odd_step_.coeffRef(18) += hip_knee_swing_angle;

  // swing hip leg
  q_even_step_.coeffRef(17) += hip_thigh_swing_angle;
  q_even_step_.coeffRef(18) += hip_knee_swing_angle;
  q_odd_step_.coeffRef(14)  += hip_thigh_swing_angle;
  q_odd_step_.coeffRef(15)  += hip_knee_swing_angle;
  step_length_ = step_length;
  v_com_ = step_length_ / t_period_;
  v_ref_.setZero();
  v_ref_.coeffRef(0) = v_com_;
  std::cout << "q_standing = " << q_standing_.transpose() << std::endl;
  std::cout << "q_even_step = " << q_even_step_.transpose() << std::endl;
  std::cout << "q_odd_step = " << q_odd_step_.transpose() << std::endl;
}


void TrottingConfigurationCost::set_q_weight(const Eigen::VectorXd& q_weight) {
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


void TrottingConfigurationCost::set_v_weight(const Eigen::VectorXd& v_weight) {
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


void TrottingConfigurationCost::set_a_weight(const Eigen::VectorXd& a_weight) {
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


void TrottingConfigurationCost::set_qf_weight(
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


void TrottingConfigurationCost::set_vf_weight(
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


double TrottingConfigurationCost::l(Robot& robot, CostFunctionData& data, 
                                    const double t, const double dtau, 
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
  return 0.5 * dtau * l;
}


double TrottingConfigurationCost::phi(Robot& robot, CostFunctionData& data, 
                                      const double t, 
                                      const SplitSolution& s) const {
  double phi = 0;
  update_q_ref(t, data.q_ref);
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, data.q_ref, data.qdiff);
    phi += (qf_weight_.array()*data.qdiff.array()*data.qdiff.array()).sum();
  }
  else {
    phi += (qf_weight_.array()*(s.q-data.q_ref).array()*(s.q-data.q_ref).array()).sum();
  }
  phi += (vf_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  return 0.5 * phi;
}


void TrottingConfigurationCost::lq(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const {
  update_q_ref(t, data.q_ref);
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


void TrottingConfigurationCost::lv(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += dtau * v_weight_.array() * (s.v.array()-v_ref_.array());
}


void TrottingConfigurationCost::la(Robot& robot, CostFunctionData& data, 
                                   const double t, const double dtau, 
                                   const SplitSolution& s, 
                                   SplitKKTResidual& kkt_residual) const {
  kkt_residual.la.array() += dtau * a_weight_.array() * s.a.array();
}


void TrottingConfigurationCost::lqq(Robot& robot, CostFunctionData& data, 
                                    const double t, const double dtau, 
                                    const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    update_q_ref(t, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += dtau * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += dtau * q_weight_;
  }
}


void TrottingConfigurationCost::lvv(Robot& robot, CostFunctionData& data, 
                                    const double t, const double dtau, 
                                    const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += dtau * v_weight_;
}


void TrottingConfigurationCost::laa(Robot& robot, CostFunctionData& data, 
                                    const double t, const double dtau, 
                                    const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qaa().diagonal().noalias() += dtau * a_weight_;
}


void TrottingConfigurationCost::phiq(Robot& robot, CostFunctionData& data, 
                                     const double t, const SplitSolution& s, 
                                     SplitKKTResidual& kkt_residual) const {
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
}


void TrottingConfigurationCost::phiv(Robot& robot, CostFunctionData& data, 
                                     const double t, const SplitSolution& s, 
                                     SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += vf_weight_.array() * (s.v.array()-v_ref_.array());
}


void TrottingConfigurationCost::phiqq(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    update_q_ref(t, data.q_ref);
    robot.dSubtractdConfigurationPlus(s.q, data.q_ref, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
  }
}


void TrottingConfigurationCost::phivv(Robot& robot, CostFunctionData& data, 
                                      const double t, const SplitSolution& s, 
                                      SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += vf_weight_;
}

} // namespace idocp