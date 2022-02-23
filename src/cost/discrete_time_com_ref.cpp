#include "robotoc/cost/discrete_time_com_ref.hpp"


namespace robotoc {

DiscreteTimeCoMRef::DiscreteTimeCoMRef(
    const std::vector<Eigen::Vector3d>& com_position_to_foot_position)
  : com_position_(), 
    com_position_to_foot_position_(com_position_to_foot_position),
    initial_com_position_(Eigen::Vector3d::Zero()),
    com_avg_(Eigen::Vector3d::Zero()),
    initial_rate_from_com_position_(0),
    has_active_contacts_(),
    has_inactive_contacts_() {
}


DiscreteTimeCoMRef::~DiscreteTimeCoMRef() {
}


void DiscreteTimeCoMRef::setCoMRef(
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const Eigen::Vector3d& initial_com_position) {
  com_position_.clear();
  has_active_contacts_.clear();
  has_inactive_contacts_.clear();
  for (int phase=0; phase<contact_sequence->numContactPhases(); ++phase) {
    const auto& contact_status = contact_sequence->contactStatus(phase);
    int num_active_contacts = 0;
    com_avg_.setZero();
    for (int i=0; i<contact_status.maxNumContacts(); ++i) {
      if (contact_status.isContactActive(i)) {
        com_avg_.noalias() += contact_status.contactPosition(i);
        com_avg_.noalias() -= com_position_to_foot_position_[i];
        ++num_active_contacts;
      }
    }
    if (num_active_contacts > 0) {
      com_avg_.array() *= (1.0/static_cast<double>(num_active_contacts));
    }
    com_position_.push_back(com_avg_);
    has_active_contacts_.push_back(contact_status.hasActiveContacts());
    has_inactive_contacts_.push_back(
        (num_active_contacts < contact_status.maxNumContacts()));
  }
  initial_com_position_ = initial_com_position;
}


void DiscreteTimeCoMRef::update_com_ref(const GridInfo& grid_info,
                                        Eigen::VectorXd& com_ref) const {
  if (has_inactive_contacts_[grid_info.contact_phase]) {
    const int next_active_contact_phase = grid_info.contact_phase + 1;
    const double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                          / static_cast<double>(grid_info.N_phase);
    com_ref = (1.0-rate) * com_position_[grid_info.contact_phase]
                + rate * com_position_[next_active_contact_phase];
  }
  else {
    com_ref = com_position_[grid_info.contact_phase];
  }
}


bool DiscreteTimeCoMRef::isActive(const GridInfo& grid_info) const {
  return true;
}

} // namespace robotoc