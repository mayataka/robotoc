#ifndef IDOCP_CONTACT_SEQUENCE_HXX
#define IDOCP_CONTACT_SEQUENCE_HXX_

#include "idocp/ocp/contact_sequence.hpp"

#include <assert.h>

namespace idocp {

inline ContactSequence::ContactSequence(const Robot& robot, const int N)
  : max_point_contacts_(robot.max_point_contacts()),
    N_(N),
    contact_sequence_(N, ContactStatus(robot.max_point_contacts())) {
}


inline ContactSequence::ContactSequence()
  : max_point_contacts_(0),
    N_(0),
    contact_sequence_() {
}


inline ContactSequence::~ContactSequence() {
}


inline void ContactSequence::activateContact(const int contact_index, 
                                             const int time_stage_begin, 
                                             const int time_stage_end) {
  assert(contact_index >= 0);
  assert(contact_index < max_point_contacts_);
  assert(time_stage_begin >= 0);
  assert(time_stage_begin < N_-1);
  assert(time_stage_end > time_stage_begin);
  assert(time_stage_end <= N_);
  for (int i=time_stage_begin; i<time_stage_end; ++i) {
    contact_sequence_[i].activateContact(contact_index);
  }
}

inline void ContactSequence::deactivateContact(const int contact_index, 
                                               const int time_stage_begin, 
                                               const int time_stage_end) {
  assert(contact_index >= 0);
  assert(contact_index < max_point_contacts_);
  assert(time_stage_begin >= 0);
  assert(time_stage_begin < N_-1);
  assert(time_stage_end > time_stage_begin);
  assert(time_stage_end <= N_);
  for (int i=time_stage_begin; i<time_stage_end; ++i) {
    contact_sequence_[i].deactivateContact(contact_index);
  }
}


inline void ContactSequence::activateContacts(
    const std::vector<int>& contact_indices, const int time_stage_begin, 
    const int time_stage_end) {
  assert(contact_indices.size() <= max_point_contacts_);
  assert(time_stage_begin >= 0);
  assert(time_stage_begin < N_);
  assert(time_stage_end > time_stage_begin);
  assert(time_stage_end <= N_);
  for (int i=time_stage_begin; i<time_stage_end; ++i) {
    contact_sequence_[i].activateContacts(contact_indices);
  }
}


inline void ContactSequence::deactivateContacts(
    const std::vector<int>& contact_indices, const int time_stage_begin, 
    const int time_stage_end) {
  assert(contact_indices.size() <= max_point_contacts_);
  assert(time_stage_begin >= 0);
  assert(time_stage_begin < N_);
  assert(time_stage_end > time_stage_begin);
  assert(time_stage_end <= N_);
  for (int i=time_stage_begin; i<time_stage_end; ++i) {
    contact_sequence_[i].deactivateContacts(contact_indices);
  }
}


inline const ContactStatus& ContactSequence::contactStatus(
    const int time_stage) const {
  assert(time_stage >= 0);
  assert(time_stage < N_);
  return contact_sequence_[time_stage];
}

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_HXX_ 