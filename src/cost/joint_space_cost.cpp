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
    vf_weight_ = qf_weight;
  }
  else {
    std::cout << "invalid argment in set_vf_weight(): size of vf_weight must be " 
              << dimv_ << std::endl;
  }
}


} // namespace idocp