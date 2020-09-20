#include "idocp/robot/point_contact.hpp"

#include <stdexcept>


namespace idocp {

PointContact::PointContact(const pinocchio::Model& model, 
                           const int contact_frame_id, 
                           const double friction_coefficient, 
                           const double restitution_coefficient) 
  : is_active_(false),
    contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    friction_coefficient_(friction_coefficient),
    restitution_coefficient_(restitution_coefficient),
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
    if (friction_coefficient_ <= 0) {
      throw std::out_of_range(
          "invalid argument: friction coefficient must be positive!");
    }
    if (restitution_coefficient < 0) {
      throw std::out_of_range(
          "invalid argument: coefficient of restitution must be nonnegative!");
    }
    if (restitution_coefficient > 1) {
      throw std::out_of_range(
          "invalid argument: coefficient of restitution must not be more than 1!");
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
    friction_coefficient_(0),
    restitution_coefficient_(0),
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


void PointContact::setFrictionCoefficient(const double friction_coefficient) {
  try {
    if (friction_coefficient <= 0) {
      throw std::out_of_range(
          "invalid argument: friction coefficient must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  friction_coefficient_ = friction_coefficient;
}


void PointContact::setRestitutionCoefficient(
    const double restitution_coefficient) {
  try {
    if (restitution_coefficient < 0) {
      throw std::out_of_range(
          "invalid argument: coefficient of restitution must be nonnegative!");
    }
    if (restitution_coefficient > 1) {
      throw std::out_of_range(
          "invalid argument: coefficient of restitution must not be more than 1!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  restitution_coefficient_ = restitution_coefficient;
}

} // namespace idocp