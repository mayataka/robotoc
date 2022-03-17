#include "robotoc/cost/discrete_time_com_ref.hpp"

#include <stdexcept>


namespace robotoc {

DiscreteTimeCoMRef::DiscreteTimeCoMRef(
    const std::vector<Eigen::Vector3d>& com_to_contact_position)
  : com_position_(), 
    com_to_contact_position_(com_to_contact_position),
    has_inactive_contacts_(),
    num_contact_phases_(1),
    first_rate_(0.0),
    last_rate_(0.0) {
}


DiscreteTimeCoMRef::~DiscreteTimeCoMRef() {
}


void DiscreteTimeCoMRef::setCoMRef(
    const std::shared_ptr<ContactSequence>& contact_sequence) {
  com_position_.clear();
  has_inactive_contacts_.clear();
  num_contact_phases_ = contact_sequence->numContactPhases();
  Eigen::Vector3d com_avg = Eigen::Vector3d::Zero();
  bool has_active_contacts_prev = true;
  for (int phase=0; phase<contact_sequence->numContactPhases(); ++phase) {
    const auto& contact_status = contact_sequence->contactStatus(phase);
    int num_active_contacts = 0;
    com_avg.setZero();
    for (int i=0; i<contact_status.maxNumContacts(); ++i) {
      if (contact_status.isContactActive(i)) {
        com_avg.noalias() += contact_status.contactPosition(i);
        com_avg.noalias() -= com_to_contact_position_[i];
        ++num_active_contacts;
      }
    }
    if (num_active_contacts > 0) {
      com_avg.array() *= (1.0/static_cast<double>(num_active_contacts));
    }
    com_position_.push_back(com_avg);
    has_inactive_contacts_.push_back(num_active_contacts < contact_status.maxNumContacts());
    if (!has_active_contacts_prev && (phase > 1)) {
      com_position_[phase-1] 
          = 0.5 * (com_position_[phase-2] + com_position_[phase]);
    }
    has_active_contacts_prev = (num_active_contacts > 0);
  }
  com_position_.push_back(com_position_.back());
}



void DiscreteTimeCoMRef::setCoMRef(
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::Vector3d& first_com_ref, const Eigen::Vector3d& last_com_ref,
    const double first_rate, const double last_rate) {
  setCoMRef(contact_sequence);
  com_position_[0] = first_com_ref;
  com_position_[contact_sequence->numContactPhases()] = last_com_ref;
  if (contact_sequence->numContactPhases() > 1) {
    if (!contact_sequence->contactStatus(1).hasActiveContacts()) {
      com_position_[1] = 0.5 * (com_position_[0] + com_position_[2]);
    }
    if (!contact_sequence->contactStatus(contact_sequence->numContactPhases()-1).hasActiveContacts()) {
      com_position_[contact_sequence->numContactPhases()-1] 
          = 0.5 * (com_position_[contact_sequence->numContactPhases()-2] 
                    + com_position_[contact_sequence->numContactPhases()]);
    }
  }
  first_rate_ = first_rate;
  last_rate_  = last_rate;
}


void DiscreteTimeCoMRef::update_com_ref(const GridInfo& grid_info,
                                        Eigen::VectorXd& com_ref) const {
  if (has_inactive_contacts_[grid_info.contact_phase]) {
    double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                    / static_cast<double>(grid_info.N_phase);
    if (grid_info.contact_phase == 0) {
      rate = first_rate_ * (1.0-rate) + rate;
    }
    else if (grid_info.contact_phase == num_contact_phases_-1) {
      rate = last_rate_ * (1.0-rate) + rate;
    }
    com_ref = (1.0-rate) * com_position_[grid_info.contact_phase]
                + rate * com_position_[grid_info.contact_phase+1];
  }
  else {
    com_ref = com_position_[grid_info.contact_phase];
  }
}


bool DiscreteTimeCoMRef::isActive(const GridInfo& grid_info) const {
  return true;
}

} // namespace robotoc