#include "robotoc/robot/surface_contact.hpp"

#include <stdexcept>


namespace robotoc {

SurfaceContact::SurfaceContact(const pinocchio::Model& model, 
                               const int contact_frame_id,
                               const double baumgarte_weight_on_velocity,
                               const double baumgarte_weight_on_position)
  : contact_frame_id_(contact_frame_id),
    parent_joint_id_(model.frames[contact_frame_id_].parent), 
    dimv_(model.nv),
    baumgarte_weight_on_velocity_(baumgarte_weight_on_velocity),
    baumgarte_weight_on_position_(baumgarte_weight_on_position), 
    jXf_(model.frames[contact_frame_id_].placement),
    X_diff_(SE3::Identity()),
    J_frame_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_v_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dq_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_dv_(Eigen::MatrixXd::Zero(6, model.nv)),
    frame_a_partial_da_(Eigen::MatrixXd::Zero(6, model.nv)),
    Jlog6_(Matrix66d::Zero()) {
  if (contact_frame_id_ < 0) {
    throw std::out_of_range(
        "[SurfaceContact] invalid argument: contact_frame_id must be non-negative!");
  }
  if (baumgarte_weight_on_velocity < 0) {
    throw std::out_of_range(
        "[SurfaceContact] invalid argument: baumgarte_weight_on_velocity must be non-negative!");
  }
  if (baumgarte_weight_on_position < 0) {
    throw std::out_of_range(
        "[SurfaceContact] invalid argument: baumgarte_weight_on_position must be non-negative!");
  }
}


SurfaceContact::SurfaceContact() 
  : contact_frame_id_(0),
    parent_joint_id_(0), 
    dimv_(0),
    baumgarte_weight_on_velocity_(0),
    baumgarte_weight_on_position_(0), 
    jXf_(),
    X_diff_(),
    J_frame_(),
    frame_v_partial_dq_(),
    frame_a_partial_dq_(),
    frame_a_partial_dv_(),
    frame_a_partial_da_(),
    Jlog6_() {
}


SurfaceContact::~SurfaceContact() {
}


void SurfaceContact::setBaumgarteWeights(
    const double baumgarte_weight_on_velocity,
    const double baumgarte_weight_on_position) {
  if (baumgarte_weight_on_velocity < 0) {
    throw std::out_of_range(
        "[SurfaceContact] invalid argument: baumgarte_weight_on_velocity must be non-negative!");
  }
  if (baumgarte_weight_on_position < 0) {
    throw std::out_of_range(
        "[SurfaceContact] invalid argument: baumgarte_weight_on_position must be non-negative!");
  }
  baumgarte_weight_on_velocity_ = baumgarte_weight_on_velocity;
  baumgarte_weight_on_position_ = baumgarte_weight_on_position;
}


void SurfaceContact::disp(std::ostream& os) const {
  os << "surface contact:" << std::endl;
  os << "  contact frame id: " << contact_frame_id_ << std::endl;
  os << "  parent joint id: " << parent_joint_id_ << std::endl;
  os << "  Baumgarte's weights on (velocity, position): (" 
     << baumgarte_weight_on_velocity_ << ", " 
     << baumgarte_weight_on_position_ << ")" << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const SurfaceContact& surface_contact) {
  surface_contact.disp(os);
  return os;
}

} // namespace robotoc 