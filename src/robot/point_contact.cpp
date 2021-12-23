#include "robotoc/robot/point_contact.hpp"

#include <stdexcept>


namespace robotoc {

PointContact::PointContact(const pinocchio::Model& model, 
                           const int contact_frame_id,
                           const double baumgarte_weight_on_velocity,
                           const double baumgarte_weight_on_position)
  : contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    baumgarte_weight_on_velocity_(baumgarte_weight_on_velocity),
    baumgarte_weight_on_position_(baumgarte_weight_on_position), 
    jXf_(model.frames[contact_frame_id_].placement),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)) {
  try {
    if (contact_frame_id_ < 0) {
      throw std::out_of_range(
          "Invalid argument: contact_frame_id must be non-negative!");
    }
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_velocity must be non-negative!");
    }
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_position must be non-negative!");
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
  : contact_frame_id_(0),
    parent_joint_id_(0), 
    dimv_(0),
    baumgarte_weight_on_velocity_(0),
    baumgarte_weight_on_position_(0), 
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


void PointContact::setBaumgarteWeights(
    const double baumgarte_weight_on_velocity,
    const double baumgarte_weight_on_position) {
  try {
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_velocity must be non-negative!");
    }
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range(
          "Invalid argument: baumgarte_weight_on_position must be non-negative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  baumgarte_weight_on_velocity_ = baumgarte_weight_on_velocity;
  baumgarte_weight_on_position_ = baumgarte_weight_on_position;
}


void PointContact::disp(std::ostream& os) const {
  os << "point contact:" << std::endl;
  os << "  contact frame id: " << contact_frame_id_ << std::endl;
  os << "  parent joint id: " << parent_joint_id_ << std::endl;
  os << "  Baumgarte's weights on (velocity, position): (" 
     << baumgarte_weight_on_velocity_ << ", " 
     << baumgarte_weight_on_position_ << ")" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const PointContact& point_contact) {
  point_contact.disp(os);
  return os;
}

} // namespace robotoc 