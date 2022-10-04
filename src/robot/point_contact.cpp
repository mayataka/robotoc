#include "robotoc/robot/point_contact.hpp"

#include <stdexcept>


namespace robotoc {

PointContact::PointContact(const pinocchio::Model& model, 
                           const ContactModelInfo& info)
  : info_(info),
    contact_frame_id_(0),
    parent_joint_id_(0), 
    dimv_(model.nv),
    jXf_(SE3::Identity()),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)) {
  if (!model.existFrame(info.frame)) {
    throw std::invalid_argument(
        "[PointContact] invalid argument: frame '" + info.frame + "' does not exit!");
  }
  contact_frame_id_ = model.getFrameId(info.frame);
  parent_joint_id_ = model.frames[contact_frame_id_].parent;
  jXf_ = model.frames[contact_frame_id_].placement;
  if (info.baumgarte_position_gain < 0) {
    throw std::out_of_range(
        "[PointContact] invalid argument: 'baumgarte_position_gain' must be non-negative!");
  }
  if (info.baumgarte_velocity_gain < 0) {
    throw std::out_of_range(
        "[PointContact] invalid argument: 'baumgarte_velocity_gain' must be non-negative!");
  }
  v_frame_.setZero();
  v_linear_skew_.setZero();
  v_angular_skew_.setZero();
}


PointContact::PointContact() 
  : info_(),
    contact_frame_id_(0),
    parent_joint_id_(0), 
    dimv_(0),
    jXf_(SE3::Identity()),
    J_frame_(),
    frame_v_partial_dq_(),
    frame_a_partial_dq_(),
    frame_a_partial_dv_(),
    frame_a_partial_da_() {
}


void PointContact::computeJointForceFromContactForce(
    const Eigen::Vector3d& contact_force, 
    pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const {
  joint_forces[parent_joint_id_] 
      = jXf_.act(pinocchio::Force(contact_force, Eigen::Vector3d::Zero()));
}


void PointContact::setBaumgarteGains(const double baumgarte_position_gain, 
                                     const double baumgarte_velocity_gain) {
  if (baumgarte_position_gain < 0) {
    throw std::out_of_range(
        "[PointContact] invalid argument: 'baumgarte_position_gain' must be non-negative!");
  }
  if (baumgarte_velocity_gain < 0) {
    throw std::out_of_range(
        "[PointContact] invalid argument: 'baumgarte_velocity_gain' must be non-negative!");
  }
  info_.baumgarte_position_gain = baumgarte_position_gain;
  info_.baumgarte_velocity_gain = baumgarte_velocity_gain;
}


const Eigen::Vector3d& PointContact::contactPosition(
    const pinocchio::Data& data) const {
  return data.oMf[contact_frame_id_].translation();
}


int PointContact::contactFrameId() const {
  return contact_frame_id_;
}


int PointContact::parentJointId() const {
  return parent_joint_id_;
}


const ContactModelInfo& PointContact::contactModelInfo() const {
  return info_;
}


void PointContact::disp(std::ostream& os) const {
  os << "PointContact:\n";
  os << "  contact frame: " << info_.frame << "\n";
  os << "  contact frame id: " << contact_frame_id_ << "\n";
  os << "  parent joint id: " << parent_joint_id_ << "\n";
  os << "  Baumgarte's gains on (position, velocity): (" 
     << info_.baumgarte_position_gain << ", " 
     << info_.baumgarte_velocity_gain << ")" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const PointContact& point_contact) {
  point_contact.disp(os);
  return os;
}

} // namespace robotoc 