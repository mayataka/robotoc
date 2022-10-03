#include "robotoc/planner/discrete_event.hpp"


namespace robotoc {

DiscreteEvent::DiscreteEvent(const ContactStatus& pre_contact_status, 
                             const ContactStatus& post_contact_status)
  : pre_contact_status_(pre_contact_status.contactTypes(), 
                        pre_contact_status.contactFrameNames()),
    post_contact_status_(pre_contact_status.contactTypes(), 
                         pre_contact_status.contactFrameNames()),
    impulse_status_(pre_contact_status.contactTypes(), 
                    pre_contact_status.contactFrameNames()),
    max_num_contacts_(pre_contact_status.maxNumContacts()),
    event_type_(DiscreteEventType::None),
    exist_impulse_(false), 
    exist_lift_(false) {
  setDiscreteEvent(pre_contact_status, post_contact_status);
}


DiscreteEvent::DiscreteEvent() 
  : pre_contact_status_(),
    post_contact_status_(),
    impulse_status_(),
    max_num_contacts_(0),
    event_type_(DiscreteEventType::None),
    exist_impulse_(false), 
    exist_lift_(false) {
}
 

void DiscreteEvent::setDiscreteEvent(
    const ContactStatus& pre_contact_status, 
    const ContactStatus& post_contact_status) {
  assert(pre_contact_status.maxNumContacts() == max_num_contacts_);
  assert(post_contact_status.maxNumContacts() == max_num_contacts_);
  exist_impulse_ = false;
  exist_lift_ = false;
  for (int i=0; i<max_num_contacts_; ++i) {
    if (pre_contact_status.isContactActive(i)) {
      impulse_status_.deactivateImpulse(i);
      if (!post_contact_status.isContactActive(i)) {
        exist_lift_ = true;
      }
    }
    else {
      if (post_contact_status.isContactActive(i)) {
        impulse_status_.activateImpulse(i);
        exist_impulse_ = true;
      }
      else {
        impulse_status_.deactivateImpulse(i);
      }
    }
  }
  impulse_status_.setFrictionCoefficients(pre_contact_status.frictionCoefficients());
  impulse_status_.setImpulseModeId(pre_contact_status.contactModeId());
  setContactPlacements(post_contact_status.contactPositions(),
                       post_contact_status.contactRotations());
  pre_contact_status_ = pre_contact_status;
  post_contact_status_ = post_contact_status;
  if (exist_impulse_) { event_type_ = DiscreteEventType::Impact; }
  else if (exist_lift_) { event_type_ = DiscreteEventType::Lift; }
  else { event_type_ = DiscreteEventType::None; }
}


void DiscreteEvent::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  impulse_status_.setContactPlacement(contact_index, contact_position);
}


void DiscreteEvent::setContactPlacement(
    const int contact_index, const Eigen::Vector3d& contact_position,
    const Eigen::Matrix3d& contact_rotation) {
  assert(contact_index >= 0);
  assert(contact_index < max_num_contacts_);
  impulse_status_.setContactPlacement(contact_index, contact_position, 
                                      contact_rotation);
}


void DiscreteEvent::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions) {
  impulse_status_.setContactPlacements(contact_positions);
}


void DiscreteEvent::setContactPlacements(
    const std::vector<Eigen::Vector3d>& contact_positions,
    const std::vector<Eigen::Matrix3d>& contact_rotations) {
  impulse_status_.setContactPlacements(contact_positions, contact_rotations);
}


void DiscreteEvent::setFrictionCoefficient(
    const int contact_index, const double friction_coefficient) {
  impulse_status_.setFrictionCoefficient(contact_index, friction_coefficient);
}


void DiscreteEvent::setFrictionCoefficients(
    const std::vector<double>& friction_coefficients) {
  impulse_status_.setFrictionCoefficients(friction_coefficients);
}


void DiscreteEvent::setContactPlacements(
    const aligned_vector<SE3>& contact_placements) {
  impulse_status_.setContactPlacements(contact_placements);
}


void DiscreteEvent::disp(std::ostream& os) const {
  auto discreteEventTypeToString = [](const DiscreteEventType& type) {
    switch (type)
    {
    case DiscreteEventType::Impact:
      return "Impact";
      break;
    case DiscreteEventType::Lift:
      return "Lift";
      break;
    default:
      return "None";
      break;
    }
  };
  os << "event type: " << discreteEventTypeToString(event_type_) << "\n";
  os << "max num contacts: " << max_num_contacts_ << "\n";
  os << "pre contact status:" << "\n";
  os << pre_contact_status_ << "\n";
  os << "post contact status:" << "\n";
  os << post_contact_status_ << "\n";
  os << "impulse status:" << "\n";
  os << impulse_status_ << std::flush;
}


std::ostream& operator<<(std::ostream& os, 
                         const DiscreteEvent& discrete_event) {
  discrete_event.disp(os);
  return os;
}

} // namespace robotoc 