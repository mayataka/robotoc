#include "robotoc/cost/discrete_time_swing_foot_ref.hpp"


namespace robotoc {

DiscreteTimeSwingFootRef::DiscreteTimeSwingFootRef(const int contact_index,
                                                   const double step_height)
  : contact_position_(),
    initial_contact_frame_position_(Eigen::Vector3d::Zero()),
    step_height_(step_height),
    initial_rate_from_step_length_(0),
    contact_index_(contact_index),
    is_contact_active_() {
}


DiscreteTimeSwingFootRef::~DiscreteTimeSwingFootRef() {
}


void DiscreteTimeSwingFootRef::setSwingFootRef(
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::Vector3d& initial_contact_frame_position, 
    const Eigen::Vector3d& contact_position_before_initial_time) {
  contact_position_.clear();
  is_contact_active_.clear();
  for (int phase=0; phase<contact_sequence->numContactPhases(); ++phase) {
    const auto& contact_status = contact_sequence->contactStatus(phase);
    contact_position_.push_back(contact_status.contactPosition(contact_index_));
    is_contact_active_.push_back(contact_status.isContactActive(contact_index_));
  }
  // Construct ref for initial swing foot using contact position before initial 
  // time of the horizon
  if (!contact_sequence->contactStatus(0).isContactActive(contact_index_)) {
    const int next_active_contact_phase = nextActiveContactPhase(0);
    Eigen::Vector3d dist = contact_position_[next_active_contact_phase] 
                            - contact_position_before_initial_time;
    const double step_length 
        = planarDistance(contact_position_[next_active_contact_phase],
                         contact_position_before_initial_time);
    const double current_step_length 
        = planarDistance(initial_contact_frame_position ,
                         contact_position_before_initial_time);
    initial_rate_from_step_length_ = current_step_length / step_length;
    initial_contact_frame_position_ = initial_contact_frame_position;
  }
}


void DiscreteTimeSwingFootRef::update_x3d_ref(const GridInfo& grid_info,
                                              Eigen::VectorXd& x3d_ref) const {
  if (!is_contact_active_[grid_info.contact_phase]) {
    const int next_active_contact_phase = nextActiveContactPhase(grid_info.contact_phase);
    const double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                          / static_cast<double>(grid_info.N_phase);
    if (grid_info.contact_phase == 0) {
      x3d_ref = (1.0-rate) * initial_contact_frame_position_
                  + rate * contact_position_[next_active_contact_phase];
      const double rate_height = (1.0-initial_rate_from_step_length_) * rate;
      if (rate_height < 0.5) {
        x3d_ref.coeffRef(2) += 2.0 * rate_height * step_height_;
      }
      else {
        x3d_ref.coeffRef(2) += 2.0 * (1.0-rate_height) * step_height_;
      }
    }
    else {
      x3d_ref = (1.0-rate) * contact_position_[grid_info.contact_phase-1]
                  + rate * contact_position_[next_active_contact_phase];
      if (rate < 0.5) {
        x3d_ref.coeffRef(2) += 2.0 * rate * step_height_;
      }
      else {
        x3d_ref.coeffRef(2) += 2.0 * (1.0-rate) * step_height_;
      }
    }
  }
}


bool DiscreteTimeSwingFootRef::isActive(const GridInfo& grid_info) const {
  return (!is_contact_active_[grid_info.contact_phase]);
}

} // namespace robotoc