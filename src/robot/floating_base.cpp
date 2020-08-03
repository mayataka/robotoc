#include "idocp/robot/floating_base.hpp"

#include <assert.h>


namespace idocp {

FloatingBase::FloatingBase(const pinocchio::Model& model) 
  : has_floating_base_(false),
    passive_joint_indices_(),
    dimv_(model.nv) {
  int total_dim_torque = -1; // Note that joint 0 is always universe.
  for (const auto& joint : model.joints) {
    if (joint.shortname() == "JointModelFreeFlyer") {
      passive_joint_indices_.push_back(total_dim_torque);
      passive_joint_indices_.push_back(total_dim_torque+1);
      passive_joint_indices_.push_back(total_dim_torque+2);
      passive_joint_indices_.push_back(total_dim_torque+3);
      passive_joint_indices_.push_back(total_dim_torque+4);
      passive_joint_indices_.push_back(total_dim_torque+5);
      total_dim_torque += 6;
      has_floating_base_ = true;
    }
    else {
      total_dim_torque += 1;
    }
  }
}


FloatingBase::FloatingBase() 
  : has_floating_base_(false),
    passive_joint_indices_(),
    dimv_(0) {
}


FloatingBase::~FloatingBase() {
}


void FloatingBase::setPassiveTorques(Eigen::VectorXd& torques) const {
  assert(torques.size() == dimv_);
  for (int i=0; i<passive_joint_indices_.size(); ++i) {
    torques.coeffRef(passive_joint_indices_[i]) = 0.0;
  }
}


void FloatingBase::computePassiveConstraintViolation(
    const Eigen::VectorXd& torques, Eigen::VectorXd& violation) const {
  assert(torques.size() == dimv_);
  assert(violation.size() == passive_joint_indices_.size());
  for (int i=0; i<passive_joint_indices_.size(); ++i) {
    violation.coeffRef(i) = torques.coeff(passive_joint_indices_[i]);
  }
}


int FloatingBase::dim_passive() const {
  return passive_joint_indices_.size();
}


std::vector<int> FloatingBase::passive_joint_indices() const {
  return passive_joint_indices_;
}


bool FloatingBase::has_floating_base() const {
  return has_floating_base_;
}

} // namespace idocp