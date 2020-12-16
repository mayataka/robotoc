#include "idocp/cost/joint_space_cost.hpp"

#include <iostream>
#include <stdexcept>


namespace idocp {

JointSpaceCost::JointSpaceCost(const Robot& robot)
  : CostFunctionComponentBase(),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    dimu_(robot.dimu()),
    q_ref_(Eigen::VectorXd::Zero(robot.dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    a_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    u_ref_(Eigen::VectorXd::Zero(robot.dimu())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    u_weight_(Eigen::VectorXd::Zero(robot.dimu())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    u_passive_weight_(Vector6d::Zero()) {
  if (robot.hasFloatingBase()) {
    robot.normalizeConfiguration(q_ref_);
  }
}


JointSpaceCost::JointSpaceCost()
  : CostFunctionComponentBase(),
    dimq_(0),
    dimv_(0),
    dimu_(),
    q_ref_(),
    v_ref_(),
    a_ref_(),
    u_ref_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    u_weight_(),
    qf_weight_(),
    vf_weight_(),
    u_passive_weight_(Vector6d::Zero()) {
}


JointSpaceCost::~JointSpaceCost() {
}


bool JointSpaceCost::useKinematics() const {
  return false;
}


void JointSpaceCost::set_q_ref(const Eigen::VectorXd& q_ref) {
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


void JointSpaceCost::set_v_ref(const Eigen::VectorXd& v_ref) {
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


void JointSpaceCost::set_a_ref(const Eigen::VectorXd& a_ref) {
  try {
    if (a_ref.size() != dimv_) {
      throw std::invalid_argument(
          "invalid size: a_ref.size() must be " + std::to_string(dimv_) + "!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  a_ref_ = a_ref;
}


void JointSpaceCost::set_u_ref(const Eigen::VectorXd& u_ref) {
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


void JointSpaceCost::set_q_weight(const Eigen::VectorXd& q_weight) {
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


void JointSpaceCost::set_v_weight(const Eigen::VectorXd& v_weight) {
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


void JointSpaceCost::set_a_weight(const Eigen::VectorXd& a_weight) {
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


void JointSpaceCost::set_u_weight(const Eigen::VectorXd& u_weight) {
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


void JointSpaceCost::set_u_passive_weight(const Vector6d& u_passive_weight) {
  u_passive_weight_ = u_passive_weight;
}


void JointSpaceCost::set_qf_weight(const Eigen::VectorXd& qf_weight) {
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


void JointSpaceCost::set_vf_weight(const Eigen::VectorXd& vf_weight) {
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


double JointSpaceCost::l(Robot& robot, CostFunctionData& data, const double t, 
                         const double dtau, const SplitSolution& s) const {
  double l = 0;
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    l += (q_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  else {
    l += (q_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  l += (v_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  l += (a_weight_.array()*(s.a-a_ref_).array()*(s.a-a_ref_).array()).sum();
  l += (u_weight_.array()*(s.u-u_ref_).array()*(s.u-u_ref_).array()).sum();
  if (robot.hasFloatingBase()) {
    l += (u_passive_weight_.array()*s.u_passive.array()*s.u_passive.array()).sum();
  }
  return 0.5 * dtau * l;
}


double JointSpaceCost::phi(Robot& robot, CostFunctionData& data, 
                           const double t, const SplitSolution& s) const {
double phi = 0;
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    phi += (qf_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  else {
    phi += (qf_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  phi += (vf_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  return 0.5 * phi;
}


void JointSpaceCost::lq(Robot& robot, CostFunctionData& data, const double t, 
                        const double dtau, const SplitSolution& s, 
                        SplitKKTResidual& kkt_residual) const {
  if (robot.hasFloatingBase()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_residual.lq().noalias()
        += dtau * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += dtau * q_weight_.array() * (s.q.array()-q_ref_.array());
  }
}


void JointSpaceCost::lv(Robot& robot, CostFunctionData& data, const double t, 
                        const double dtau, const SplitSolution& s, 
                        SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += dtau * v_weight_.array() * (s.v.array()-v_ref_.array());
}


void JointSpaceCost::la(Robot& robot, CostFunctionData& data, const double t, 
                        const double dtau, const SplitSolution& s, 
                        SplitKKTResidual& kkt_residual) const {
  kkt_residual.la.array()
      += dtau * a_weight_.array() * (s.a.array()-a_ref_.array());
}


void JointSpaceCost::lu(Robot& robot, CostFunctionData& data, const double t, 
                        const double dtau, const SplitSolution& s, 
                        SplitKKTResidual& kkt_residual) const {
  if (robot.hasFloatingBase()) {
    kkt_residual.lu_passive.array()
      += dtau * u_passive_weight_.array() * s.u_passive.array();
  }
  kkt_residual.lu().array() 
      += dtau * u_weight_.array() * (s.u.array()-u_ref_.array());
}


void JointSpaceCost::lqq(Robot& robot, CostFunctionData& data, const double t, 
                         const double dtau, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += dtau * data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += dtau * q_weight_;
  }
}


void JointSpaceCost::lvv(Robot& robot, CostFunctionData& data, const double t, 
                         const double dtau, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += dtau * v_weight_;
}


void JointSpaceCost::laa(Robot& robot, CostFunctionData& data, const double t, 
                         const double dtau, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qaa().diagonal().noalias() += dtau * a_weight_;
}


void JointSpaceCost::luu(Robot& robot, CostFunctionData& data, const double t, 
                         const double dtau, const SplitSolution& s, 
                         SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    kkt_matrix.Quu_passive_topLeft().diagonal().noalias()
        += dtau * u_passive_weight_;
  }
  kkt_matrix.Quu().diagonal().noalias() += dtau * u_weight_;
}


void JointSpaceCost::phiq(Robot& robot, CostFunctionData& data, const double t, 
                          const SplitSolution& s,
                          SplitKKTResidual& kkt_residual) const {
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
}


void JointSpaceCost::phiv(Robot& robot, CostFunctionData& data, const double t, 
                          const SplitSolution& s,
                          SplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += vf_weight_.array() * (s.v.array()-v_ref_.array());
}


void JointSpaceCost::phiqq(Robot& robot, CostFunctionData& data, const double t, 
                           const SplitSolution& s, 
                           SplitKKTMatrix& kkt_matrix) const {
  if (robot.hasFloatingBase()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * qf_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += qf_weight_;
  }
}


void JointSpaceCost::phivv(Robot& robot, CostFunctionData& data, const double t, 
                           const SplitSolution& s, 
                           SplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += vf_weight_;
}

} // namespace idocp