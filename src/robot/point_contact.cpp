#include "idocp/robot/point_contact.hpp"

#include <stdexcept>


namespace idocp {

PointContact::PointContact(const pinocchio::Model& model, 
                           const int contact_frame_id, 
                           const double baumgarte_weight_on_velocity, 
                           const double baumgarte_weight_on_position)
  : is_active_(false),
    contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    baumgarte_weight_on_velocity_(baumgarte_weight_on_velocity),
    baumgarte_weight_on_position_(baumgarte_weight_on_position),
    contact_point_(Eigen::Vector3d::Zero()),
    jXf_(model.frames[contact_frame_id_].placement),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)) {
  try {
    if (contact_frame_id_ < 0) {
      throw std::out_of_range("invalid argument: contct frame index must be nonnegative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  try {
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range("invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  try {
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range("invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  v_frame_.setZero();
  v_linear_skew_.setZero();
  v_angular_skew_.setZero();
}


PointContact::PointContact() 
  : is_active_(false),
    contact_frame_id_(0),
    parent_joint_id_(0), 
    dimv_(0),
    baumgarte_weight_on_velocity_(0),
    baumgarte_weight_on_position_(0),
    contact_point_(Eigen::Vector3d::Zero()),
    jXf_(),
    J_frame_(),
    frame_v_partial_dq_(),
    frame_a_partial_dq_(),
    frame_a_partial_dv_(),
    frame_a_partial_da_() {
  v_frame_.setZero();
  v_linear_skew_.setZero();
  v_angular_skew_.setZero();
}


PointContact::~PointContact() {
}


void PointContact::resetBaugrarteParameters(
    const double baumgarte_weight_on_velocity, 
    const double baumgarte_weight_on_position) {
  try {
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range("invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  try {
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range("invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  baumgarte_weight_on_velocity_ = baumgarte_weight_on_velocity;
  baumgarte_weight_on_position_ = baumgarte_weight_on_position;
}


void PointContact::resetContactPoint(const Eigen::Vector3d& contact_point) {
  contact_point_ = contact_point;
}


void PointContact::resetContactPointByCurrentKinematics(
    const pinocchio::Data& data) {
  contact_point_ = data.oMf[contact_frame_id_].translation();
}


void PointContact::computeJointForceFromContactForce(
    const Eigen::Vector3d& contact_force, 
    pinocchio::container::aligned_vector<pinocchio::Force>& joint_forces) const {
  joint_forces[parent_joint_id_] 
      = jXf_.act(pinocchio::Force(contact_force, Eigen::Vector3d::Zero()));
}


void PointContact::activate() {
  is_active_ = true;
}


void PointContact::deactivate() {
  is_active_ = false;
}


bool PointContact::isActive() const {
  return is_active_;
}


int PointContact::contact_frame_id() const {
  return contact_frame_id_;
}


int PointContact::parent_joint_id() const {
  return parent_joint_id_;
}


double PointContact::baumgarte_weight_on_velocity() const {
  return baumgarte_weight_on_velocity_;
}


double PointContact::baumgarte_weight_on_position() const {
  return baumgarte_weight_on_position_;
}


Eigen::Vector3d PointContact::contact_point() const {
  return contact_point_;
}

} // namespace idocp