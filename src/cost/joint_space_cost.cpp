#include "idocp/cost/joint_space_cost.hpp"

#include <assert.h>


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
  assert(q_ref.size() == dimq_);
  q_ref_ = q_ref;
}


void JointSpaceCost::set_v_ref(const Eigen::VectorXd& v_ref) {
  assert(v_ref.size() == dimv_);
  v_ref_ = v_ref;
}


void JointSpaceCost::set_a_ref(const Eigen::VectorXd& a_ref) {
  assert(a_ref.size() == dimv_);
  a_ref_ = a_ref;
}


void JointSpaceCost::set_u_ref(const Eigen::VectorXd& u_ref) {
  assert(u_ref.size() == dimv_);
  u_ref_ = u_ref;
}


void JointSpaceCost::set_q_weight(const Eigen::VectorXd& q_weight) {
  assert(q_weight.size() == dimv_);
  q_weight_ = q_weight;
}


void JointSpaceCost::set_v_weight(const Eigen::VectorXd& v_weight) {
  assert(v_weight.size() == dimv_);
  v_weight_ = v_weight;
}


void JointSpaceCost::set_a_weight(const Eigen::VectorXd& a_weight) {
  assert(a_weight.size() == dimv_);
  a_weight_ = a_weight;
}


void JointSpaceCost::set_u_weight(const Eigen::VectorXd& u_weight) {
  assert(u_weight.size() == dimv_);
  u_weight_ = u_weight;
}


void JointSpaceCost::set_qf_weight(const Eigen::VectorXd& qf_weight) {
  assert(qf_weight.size() == dimv_);
  qf_weight_ = qf_weight;
}


void JointSpaceCost::set_vf_weight(const Eigen::VectorXd& vf_weight) {
  assert(vf_weight.size() == dimv_);
  vf_weight_ = vf_weight;
}


double JointSpaceCost::l(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const Eigen::Ref<const Eigen::VectorXd>& q, 
                         const Eigen::Ref<const Eigen::VectorXd>& v, 
                         const Eigen::Ref<const Eigen::VectorXd>& a, 
                         const Eigen::Ref<const Eigen::VectorXd>& f, 
                         const Eigen::Ref<const Eigen::VectorXd>& u) const {
  double l = 0;
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(q, q_ref_, data.q_diff);
    l += (q_weight_.array()*(data.q_diff).array()*(data.q_diff).array()).sum();
  }
  else {
    l += (q_weight_.array()*(q-q_ref_).array()*(q-q_ref_).array()).sum();
  }
  l += (v_weight_.array()*(v-v_ref_).array()*(v-v_ref_).array()).sum();
  l += (a_weight_.array()*(a-a_ref_).array()*(a-a_ref_).array()).sum();
  l += (u_weight_.array()*(u-u_ref_).array()*(u-u_ref_).array()).sum();
  return 0.5 * dtau * l;
}


double JointSpaceCost::phi(const Robot& robot, CostFunctionData& data, 
                           const double t, 
                           const Eigen::Ref<const Eigen::VectorXd>& q, 
                           const Eigen::Ref<const Eigen::VectorXd>& v) const {
  double phi = 0;
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(q, q_ref_, data.q_diff);
    phi += (qf_weight_.array()*(data.q_diff).array()*(data.q_diff).array()).sum();
  }
  else {
    phi += (qf_weight_.array()*(q-q_ref_).array()*(q-q_ref_).array()).sum();
  }
  phi += (vf_weight_.array()*(v-v_ref_).array()*(v-v_ref_).array()).sum();
  return 0.5 * phi;
}


void JointSpaceCost::lq(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const Eigen::Ref<const Eigen::VectorXd>& q, 
                        const Eigen::Ref<const Eigen::VectorXd>& v, 
                        const Eigen::Ref<const Eigen::VectorXd>& a, 
                        Eigen::Ref<Eigen::VectorXd> lq) const {
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(q, q_ref_, data.q_diff);
    robot.dSubtractdConfigurationPlus(q, q_ref_, data.Jq_diff);
    lq = dtau * data.Jq_diff.transpose() * q_weight_.asDiagonal() * data.q_diff;
  }
  else {
    lq.array() = dtau * q_weight_.array() * (q.array()-q_ref_.array());
  }
}


void JointSpaceCost::lv(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const Eigen::Ref<const Eigen::VectorXd>& q, 
                        const Eigen::Ref<const Eigen::VectorXd>& v, 
                        const Eigen::Ref<const Eigen::VectorXd>& a, 
                        Eigen::Ref<Eigen::VectorXd> lv) const {
  lv.array() = dtau * v_weight_.array() * (v.array()-v_ref_.array());
}


void JointSpaceCost::la(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const Eigen::Ref<const Eigen::VectorXd>& q, 
                        const Eigen::Ref<const Eigen::VectorXd>& v, 
                        const Eigen::Ref<const Eigen::VectorXd>& a, 
                        Eigen::Ref<Eigen::VectorXd> la) const {
  la.array() = dtau * a_weight_.array() * (a.array()-a_ref_.array());
}


void JointSpaceCost::lu(const Robot& robot, CostFunctionData& data, 
                        const double t, const double dtau, 
                        const Eigen::Ref<const Eigen::VectorXd>& u, 
                        Eigen::Ref<Eigen::VectorXd> lu) const {
  lu.array() = dtau * u_weight_.array() * (u.array()-u_ref_.array());
}


void JointSpaceCost::lqq(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const Eigen::Ref<const Eigen::VectorXd>& q, 
                         const Eigen::Ref<const Eigen::VectorXd>& v, 
                         const Eigen::Ref<const Eigen::VectorXd>& a, 
                         Eigen::Ref<Eigen::MatrixXd> lqq) const {
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationPlus(q, q_ref_, data.Jq_diff);
    lqq = dtau * data.Jq_diff.transpose() 
               * q_weight_.asDiagonal() * data.Jq_diff;
  }
  else {
    for (int i=0; i<robot.dimv(); ++i) {
      lqq.coeffRef(i, i) = dtau * q_weight_.coeff(i);
    }
  }
}


void JointSpaceCost::lvv(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const Eigen::Ref<const Eigen::VectorXd>& q, 
                         const Eigen::Ref<const Eigen::VectorXd>& v, 
                         const Eigen::Ref<const Eigen::VectorXd>& a, 
                         Eigen::Ref<Eigen::MatrixXd> lvv) const {
  for (int i=0; i<robot.dimv(); ++i) {
    lvv.coeffRef(i, i) = dtau * v_weight_.coeff(i);
  }
}


void JointSpaceCost::laa(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const Eigen::Ref<const Eigen::VectorXd>& q, 
                         const Eigen::Ref<const Eigen::VectorXd>& v, 
                         const Eigen::Ref<const Eigen::VectorXd>& a, 
                         Eigen::Ref<Eigen::MatrixXd> laa) const {
  for (int i=0; i<robot.dimv(); ++i) {
    laa.coeffRef(i, i) = dtau * a_weight_.coeff(i);
  }
}


void JointSpaceCost::luu(const Robot& robot, CostFunctionData& data, 
                         const double t, const double dtau, 
                         const Eigen::Ref<const Eigen::VectorXd>& u, 
                         Eigen::Ref<Eigen::MatrixXd> luu) const {
  for (int i=0; i<robot.dimv(); ++i) {
    luu.coeffRef(i, i) = dtau * u_weight_.coeff(i);
  }
}


void JointSpaceCost::augment_lqq(const Robot& robot, CostFunctionData& data, 
                                 const double t, const double dtau, 
                                 const Eigen::Ref<const Eigen::VectorXd>& q, 
                                 const Eigen::Ref<const Eigen::VectorXd>& v, 
                                 const Eigen::Ref<const Eigen::VectorXd>& a, 
                                 Eigen::Ref<Eigen::MatrixXd> lqq) const {
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationPlus(q, q_ref_, data.Jq_diff);
    lqq.noalias() += dtau * data.Jq_diff.transpose() 
                          * q_weight_.asDiagonal() * data.Jq_diff;
  }
  else {
    for (int i=0; i<robot.dimv(); ++i) {
      lqq.coeffRef(i, i) += dtau * q_weight_.coeff(i);
    }
  }
}


void JointSpaceCost::augment_lvv(const Robot& robot, CostFunctionData& data, 
                                 const double t, const double dtau, 
                                 const Eigen::Ref<const Eigen::VectorXd>& q, 
                                 const Eigen::Ref<const Eigen::VectorXd>& v, 
                                 const Eigen::Ref<const Eigen::VectorXd>& a, 
                                 Eigen::Ref<Eigen::MatrixXd> lvv) const {
  for (int i=0; i<robot.dimv(); ++i) {
    lvv.coeffRef(i, i) += dtau * v_weight_.coeff(i);
  }
}


void JointSpaceCost::augment_laa(const Robot& robot, CostFunctionData& data, 
                                 const double t, const double dtau, 
                                 const Eigen::Ref<const Eigen::VectorXd>& q, 
                                 const Eigen::Ref<const Eigen::VectorXd>& v, 
                                 const Eigen::Ref<const Eigen::VectorXd>& a, 
                                 Eigen::Ref<Eigen::MatrixXd> laa) const {
  for (int i=0; i<robot.dimv(); ++i) {
    laa.coeffRef(i, i) += dtau * a_weight_.coeff(i);
  }
}


void JointSpaceCost::augment_luu(const Robot& robot, CostFunctionData& data, 
                                 const double t, const double dtau, 
                                 const Eigen::Ref<const Eigen::VectorXd>& u, 
                                 Eigen::Ref<Eigen::MatrixXd> luu) const {
  for (int i=0; i<robot.dimv(); ++i) {
    luu.coeffRef(i, i) += dtau * u_weight_.coeff(i);
  }
}


void JointSpaceCost::phiq(const Robot& robot, CostFunctionData& data, 
                          const double t, 
                          const Eigen::Ref<const Eigen::VectorXd>& q, 
                          const Eigen::Ref<const Eigen::VectorXd>& v, 
                          Eigen::Ref<Eigen::VectorXd> phiq) const {
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(q, q_ref_, data.q_diff);
    robot.dSubtractdConfigurationPlus(q, q_ref_, data.Jq_diff);
    phiq = data.Jq_diff.transpose() * qf_weight_.asDiagonal() * data.q_diff;
  }
  else {
    phiq.array() = qf_weight_.array() * (q.array()-q_ref_.array());
  }
}


void JointSpaceCost::phiv(const Robot& robot, CostFunctionData& data, 
                          const double t, 
                          const Eigen::Ref<const Eigen::VectorXd>& q, 
                          const Eigen::Ref<const Eigen::VectorXd>& v, 
                          Eigen::Ref<Eigen::VectorXd> phiv) const {
  phiv.array() = vf_weight_.array() * (v.array()-v_ref_.array());
}


void JointSpaceCost::phiqq(const Robot& robot, CostFunctionData& data, 
                           const double t, 
                           const Eigen::Ref<const Eigen::VectorXd>& q, 
                           const Eigen::Ref<const Eigen::VectorXd>& v, 
                           Eigen::Ref<Eigen::MatrixXd> phiqq) const {
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationPlus(q, q_ref_, data.Jq_diff);
    phiqq = data.Jq_diff.transpose() * qf_weight_.asDiagonal() * data.Jq_diff;
  }
  else {
    for (int i=0; i<robot.dimv(); ++i) {
      phiqq.coeffRef(i, i) = qf_weight_.coeff(i);
    }
  }
}


void JointSpaceCost::phivv(const Robot& robot, CostFunctionData& data, 
                           const double t, 
                           const Eigen::Ref<const Eigen::VectorXd>& q, 
                           const Eigen::Ref<const Eigen::VectorXd>& v, 
                           Eigen::Ref<Eigen::MatrixXd> phivv) const {
  for (int i=0; i<robot.dimv(); ++i) {
    phivv.coeffRef(i, i) = vf_weight_.coeff(i);
  }
}

} // namespace idocp