#include "idocp/cost/joint_space_impulse_cost.hpp"

#include <iostream>


namespace idocp {

JointSpaceImpulseCost::JointSpaceImpulseCost(const Robot& robot)
  : ImpulseCostFunctionComponentBase(),
    dimq_(robot.dimq()),
    dimv_(robot.dimv()),
    q_ref_(Eigen::VectorXd::Zero(robot.dimq())),
    v_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    dv_ref_(Eigen::VectorXd::Zero(robot.dimv())),
    q_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    v_weight_(Eigen::VectorXd::Zero(robot.dimv())),
    dv_weight_(Eigen::VectorXd::Zero(robot.dimv())) {
  if (robot.has_floating_base()) {
    robot.normalizeConfiguration(q_ref_);
  }
}


JointSpaceImpulseCost::JointSpaceImpulseCost()
  : ImpulseCostFunctionComponentBase(),
    dimq_(0),
    dimv_(0),
    q_ref_(),
    v_ref_(),
    dv_ref_(),
    q_weight_(),
    v_weight_(),
    dv_weight_() {
}


JointSpaceImpulseCost::~JointSpaceImpulseCost() {
}


void JointSpaceImpulseCost::set_q_ref(const Eigen::VectorXd& q_ref) {
  if (q_ref.size() == dimq_) {
    q_ref_ = q_ref;
  }
  else {
    std::cout << "invalid argment in set_q_ref(): size of q_ref must be " 
              << dimq_ << std::endl;
  }
}


void JointSpaceImpulseCost::set_v_ref(const Eigen::VectorXd& v_ref) {
  if (v_ref.size() == dimv_) {
    v_ref_ = v_ref;
  }
  else {
    std::cout << "invalid argment in set_v_ref(): size of v_ref must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceImpulseCost::set_dv_ref(const Eigen::VectorXd& dv_ref) {
  if (dv_ref.size() == dimv_) {
    dv_ref_ = dv_ref;
  }
  else {
    std::cout << "invalid argment in set_v_ref(): size of v_ref must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceImpulseCost::set_q_weight(const Eigen::VectorXd& q_weight) {
  if (q_weight.size() == dimv_) {
    q_weight_ = q_weight;
  }
  else {
    std::cout << "invalid argment in set_q_weight(): size of q_weight must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceImpulseCost::set_v_weight(const Eigen::VectorXd& v_weight) {
  if (v_weight.size() == dimv_) {
    v_weight_ = v_weight;
  }
  else {
    std::cout << "invalid argment in set_v_weight(): size of v_weight must be " 
              << dimv_ << std::endl;
  }
}


void JointSpaceImpulseCost::set_dv_weight(const Eigen::VectorXd& dv_weight) {
  if (dv_weight.size() == dimv_) {
    dv_weight_ = dv_weight;
  }
  else {
    std::cout << "invalid argment in set_v_weight(): size of v_weight must be " 
              << dimv_ << std::endl;
  }
}


double JointSpaceImpulseCost::l(Robot& robot, CostFunctionData& data, 
                                const double t, 
                                const ImpulseSplitSolution& s) const {
  double l = 0;
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    l += (q_weight_.array()*(data.qdiff).array()*(data.qdiff).array()).sum();
  }
  else {
    l += (q_weight_.array()*(s.q-q_ref_).array()*(s.q-q_ref_).array()).sum();
  }
  l += (v_weight_.array()*(s.v-v_ref_).array()*(s.v-v_ref_).array()).sum();
  l += (dv_weight_.array()*(s.dv-dv_ref_).array()*(s.dv-dv_ref_).array()).sum();
  return 0.5 * l;
}


void JointSpaceImpulseCost::lq(Robot& robot, CostFunctionData& data, 
                               const double t, const ImpulseSplitSolution& s, 
                               ImpulseSplitKKTResidual& kkt_residual) const {
  if (robot.has_floating_base()) {
    robot.subtractConfiguration(s.q, q_ref_, data.qdiff);
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_residual.lq().noalias()
        += data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.qdiff;
  }
  else {
    kkt_residual.lq().array()
        += q_weight_.array() * (s.q.array()-q_ref_.array());
  }
}


void JointSpaceImpulseCost::lv(Robot& robot, CostFunctionData& data, 
                               const double t, const ImpulseSplitSolution& s, 
                               ImpulseSplitKKTResidual& kkt_residual) const {
  kkt_residual.lv().array()
      += v_weight_.array() * (s.v.array()-v_ref_.array());
}



void JointSpaceImpulseCost::ldv(Robot& robot, CostFunctionData& data, 
                                const double t, const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTResidual& kkt_residual) const {
  kkt_residual.ldv.array() 
      += dv_weight_.array() * (s.dv.array()-dv_ref_.array());
}



void JointSpaceImpulseCost::lqq(Robot& robot, CostFunctionData& data, 
                                const double t, const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTMatrix& kkt_matrix) const {
  if (robot.has_floating_base()) {
    robot.dSubtractdConfigurationPlus(s.q, q_ref_, data.J_qdiff);
    kkt_matrix.Qqq().noalias()
        += data.J_qdiff.transpose() * q_weight_.asDiagonal() * data.J_qdiff;
  }
  else {
    kkt_matrix.Qqq().diagonal().noalias() += q_weight_;
  }
}


void JointSpaceImpulseCost::lvv(Robot& robot, CostFunctionData& data, 
                                const double t, const ImpulseSplitSolution& s, 
                                ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qvv().diagonal().noalias() += v_weight_;
}


void JointSpaceImpulseCost::ldvdv(Robot& robot, CostFunctionData& data, 
                                  const double t, const ImpulseSplitSolution& s, 
                                  ImpulseSplitKKTMatrix& kkt_matrix) const {
  kkt_matrix.Qdvdv().diagonal().noalias() += dv_weight_;
}

} // namespace idocp