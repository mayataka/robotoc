#include "idocp/robot/point_contact.hpp"

#include <stdexcept>


namespace idocp {

PointContact::PointContact(const pinocchio::Model& model, 
                           const int contact_frame_id, 
                           const double mu,
                           const double baumgarte_weight_on_velocity, 
                           const double baumgarte_weight_on_position)
  : is_active_(false),
    contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    mu_(mu),
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
      throw std::out_of_range(
          "invalid argument: contct frame index must be nonnegative!");
    }
    if (mu <= 0) {
      throw std::out_of_range(
          "invalid argument: friction coefficient mu must be positive!");
    }
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range(
          "invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
    }
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range(
          "invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
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
    mu_(0),
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


void PointContact::setFrictionCoefficient(const double mu) {
  try {
    if (mu <= 0) {
      throw std::out_of_range(
          "invalid argument: friction coefficient mu must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  mu_ = mu;
}


void PointContact::setBaugrarteParameters(
    const double baumgarte_weight_on_velocity, 
    const double baumgarte_weight_on_position) {
  try {
    if (baumgarte_weight_on_velocity < 0) {
      throw std::out_of_range(
          "invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  try {
    if (baumgarte_weight_on_position < 0) {
      throw std::out_of_range(
          "invalid argument: weight on Baumgarte's stabilization must be nonnegative!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  baumgarte_weight_on_velocity_ = baumgarte_weight_on_velocity;
  baumgarte_weight_on_position_ = baumgarte_weight_on_position;
}

} // namespace idocp