#include "idocp/cost/joint_space_cost.hpp"

#include <iostream>


namespace idocp {

JointSpaceCost::JointSpaceCost(const Robot& robot)
  : CostFunctionComponentBase(),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    q_ref_(Eigen::VectorXd::Zero(robot.dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    a_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    u_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    a_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    u_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    qf_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    vf_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
}


JointSpaceCost::JointSpaceCost()
  : CostFunctionComponentBase(),
    dimq_(0),
    dimv_(0),
    q_ref_(),
    v_ref_(),
    a_ref_(),
    u_ref_(),
    q_weight_(),
    v_weight_(),
    a_weight_(),
    u_weight_(),
    qf_weight_(),
    vf_weight_() {
}


JointSpaceCost::~JointSpaceCost() {
}


void JointSpaceCost::set_q_ref(const Eigen::VectorXd& q_ref) {
  if (q_ref.size() == dimq_) {
    q_ref_ = q_ref;
  }
  else {
    std::cout << "invalid argment in set_q_ref(): size of q_ref must be " 
              << dimq_ << std::endl;
  }
}


void JointSpaceCost::set_v_ref(const Eigen::VectorXd& v_ref) {
  if (v_ref.size() == dimv_) {
    v_ref_ = v_ref;
  }
  else {
    std::cout << "invalid argment in set_v_ref(): size of v_ref must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_a_ref(const Eigen::VectorXd& a_ref) {
  if (a_ref.size() == dimv_) {
    a_ref_ = a_ref;
  }
  else {
    std::cout << "invalid argment in set_a_ref(): size of a_ref must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_u_ref(const Eigen::VectorXd& u_ref) {
  if (u_ref.size() == dimv_) {
    u_ref_ = u_ref;
  }
  else {
    std::cout << "invalid argment in set_u_ref(): size of u_ref must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_q_weight(const Eigen::VectorXd& q_weight) {
  if (q_weight.size() == dimv_) {
    q_weight_ = q_weight;
  }
  else {
    std::cout << "invalid argment in set_q_weight(): size of q_weight must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_v_weight(const Eigen::VectorXd& v_weight) {
  if (v_weight.size() == dimv_) {
    v_weight_ = v_weight;
  }
  else {
    std::cout << "invalid argment in set_v_weight(): size of v_weight must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_a_weight(const Eigen::VectorXd& a_weight) {
  if (a_weight.size() == dimv_) {
    a_weight_ = a_weight;
  }
  else {
    std::cout << "invalid argment in set_a_weight(): size of a_weight must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_u_weight(const Eigen::VectorXd& u_weight) {
  if (u_weight.size() == dimv_) {
    u_weight_ = u_weight;
  }
  else {
    std::cout << "invalid argment in set_u_weight(): size of u_weight must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_qf_weight(const Eigen::VectorXd& qf_weight) {
  if (qf_weight.size() == dimv_) {
    qf_weight_ = qf_weight;
  }
  else {
    std::cout << "invalid argment in set_qf_weight(): size of qf_weight must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceCost::set_vf_weight(const Eigen::VectorXd& vf_weight) {
  if (vf_weight.size() == dimv_) {
    vf_weight_ = vf_weight;
  }
  else {
    std::cout << "invalid argment in set_vf_weight(): size of vf_weight must be " 
              << dimv_ << std::endl;
  }
}


double JointSpaceCost::l(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s) const {
  double l = 0;
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
    l += (q_weight_.array()*(data.q_diff).array()*(data.q_diff).array()).sum();
  }
  else {
    l += (q_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  l += (v_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  l += (a_weight_.array()*(s.a-a_ref_).array()*(s.a-a_ref_).array()).sum();
  l += (u_weight_.array()*(s.u-u_ref_).array()*(s.u-u_ref_).array()).sum();
  return 0.5 * dtau * l;
}


double JointSpaceCost::phi(const Robot& robot, CostFunctionData& data, 
                           const double t, const SplitSolution& s) const {
  double phi = 0;
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
    phi += (qf_weight_.array()*(data.q_diff).array()*(data.q_diff).array()).sum();
  }
  else {
    phi += (qf_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  phi += (vf_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  return 0.5 * phi;
}


void JointSpaceCost::lq(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const SplitSolution& s, 
                        KKTResidual& kkt_residual) const {
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
    kkt_residual.lq().noalias()
        += dtau * data.Jq_diff.transpose() * q_weight_.asDiagonal() * data.q_diff;
  }
  else {
    kkt_residual.lq().array()
        += dtau * q_weight_.array() * (s.q.array()-q_ref_.array());
  }
}


void JointSpaceCost::lv(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const SplitSolution& s, 
                        KKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += dtau * v_weight_.array() * (s.v.array()-v_ref_.array());
}


void JointSpaceCost::la(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const SplitSolution& s, 
                        KKTResidual& kkt_residual) const {
  kkt_residual.la().array()
      += dtau * a_weight_.array() * (s.a.array()-a_ref_.array());
}


void JointSpaceCost::lf(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const SplitSolution& s, 
                        KKTResidual& kkt_residual) const {
  // do nothing
}


void JointSpaceCost::lu(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const SplitSolution& s, 
                        KKTResidual& kkt_residual) const {
  kkt_residual.lu.array()
      += dtau * u_weight_.array() * (s.u.array()-u_ref_.array());
}


void JointSpaceCost::lqq(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
    kkt_matrix.Qqq().noalias()
        += dtau * data.Jq_diff.transpose() * q_weight_.asDiagonal() * data.Jq_diff;
  }
  else {
    kkt_matrix.Qqq() += dtau * q_weight_.asDiagonal();
  }
}


void JointSpaceCost::lvv(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv() += dtau * v_weight_.asDiagonal();
}


void JointSpaceCost::laa(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  kkt_matrix.Qaa() += dtau * a_weight_.asDiagonal();
}


void JointSpaceCost::lff(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  // do nothing
}


void JointSpaceCost::luu(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const SplitSolution& s, KKTMatrix& kkt_matrix) const {
  kkt_matrix.luu += dtau * u_weight_.asDiagonal();
}


void JointSpaceCost::phiq(const Robot& robot, CostFunctionData& data, 
                          const double t, const SplitSolution& s,
                          KKTResidual& kkt_residual) const {
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, q_ref_, data.q_diff);
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
    kkt_residual.phiq.noalias()
        += data.Jq_diff.transpose() * qf_weight_.asDiagonal() * data.q_diff;
  }
  else {
    kkt_residual.phiq.array()
        += qf_weight_.array() * (s.q.array()-q_ref_.array());
  }
}


void JointSpaceCost::phiv(const Robot& robot, CostFunctionData& data, 
                          const double t, const SplitSolution& s,
                          KKTResidual& kkt_residual) const {
    kkt_residual.phiv.array()
        += vf_weight_.array() * (s.v.array()-v_ref_.array());
}


void JointSpaceCost::phiqq(const Robot& robot, CostFunctionData& data, 
                           const double t, const SplitSolution& s,
                           KKTMatrix& kkt_matrix) const {
    if (robot.has_floating_base()) {
      robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.Jq_diff);
      kkt_matrix.Qqq().noalias()
          += data.Jq_diff.transpose() * qf_weight_.asDiagonal() * data.Jq_diff;
    }
    else {
      kkt_matrix.Qqq() += qf_weight_.asDiagonal();
    }
}


void JointSpaceCost::phivv(const Robot& robot, CostFunctionData& data, 
                          const double t, const SplitSolution& s,
                          KKTMatrix& kkt_matrix) const {

    kkt_matrix.Qvv() += vf_weight_.asDiagonal();
}

} // namespace idocp