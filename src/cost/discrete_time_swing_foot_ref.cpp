#include "robotoc/cost/discrete_time_swing_foot_ref.hpp"


namespace robotoc {

DiscreteTimeSwingFootRef::DiscreteTimeSwingFootRef(const int contact_index,
                                                   const double step_height)
  : contact_index_(contact_index),
    num_contact_phases_(1),
    step_height_(step_height),
    first_rate_(0.0),
    last_rate_(0.0),
    contact_position_(),
    is_contact_active_() {
}


DiscreteTimeSwingFootRef::~DiscreteTimeSwingFootRef() {
}


void DiscreteTimeSwingFootRef::setSwingFootRef(
    const std::shared_ptr<ContactSequence>& contact_sequence) {
  contact_position_.clear();
  is_contact_active_.clear();
  num_contact_phases_ = contact_sequence->numContactPhases();
  for (int phase=0; phase<contact_sequence->numContactPhases(); ++phase) {
    const auto& contact_status = contact_sequence->contactStatus(phase);
    contact_position_.push_back(contact_status.contactPosition(contact_index_));
    is_contact_active_.push_back(contact_status.isContactActive(contact_index_));
  }
  contact_position_.push_back(contact_position_.back());
  first_rate_ = 1.0;
  last_rate_ = 1.0;
}


void DiscreteTimeSwingFootRef::setSwingFootRef(
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::Vector3d& first_contact_position, 
    const Eigen::Vector3d& last_contact_position, 
    const double first_rate, const double last_rate) {
  setSwingFootRef(contact_sequence);
  contact_position_[0] = first_contact_position;
  contact_position_[contact_sequence->numContactPhases()] = last_contact_position;
  first_rate_ = first_rate;
  last_rate_  = last_rate;
}


void DiscreteTimeSwingFootRef::updateRef(const GridInfo& grid_info,
                                         Eigen::VectorXd& x3d_ref) const {
  if (!is_contact_active_[grid_info.contact_phase]) {
    double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                    / static_cast<double>(grid_info.N_phase);
    if (grid_info.contact_phase == 0) {
      rate = first_rate_ * (1.0-rate) + rate;
    }
    else if (grid_info.contact_phase == num_contact_phases_-1) {
      rate = last_rate_ * (1.0-rate) + rate;
    }
    if (grid_info.contact_phase == 0 ) {
      x3d_ref = (1.0-rate) * contact_position_[0] + rate * contact_position_[1];
    }
    else {
      x3d_ref = (1.0-rate) * contact_position_[grid_info.contact_phase-1]
                  + rate * contact_position_[grid_info.contact_phase+1];
    }
    if (rate < 0.5) {
      x3d_ref.coeffRef(2) += 2.0 * rate * step_height_;
    }
    else {
      x3d_ref.coeffRef(2) += 2.0 * (1.0-rate) * step_height_;
    }
  }
}


bool DiscreteTimeSwingFootRef::isActive(const GridInfo& grid_info) const {
  return (!is_contact_active_[grid_info.contact_phase]);
}

} // namespace robotoc