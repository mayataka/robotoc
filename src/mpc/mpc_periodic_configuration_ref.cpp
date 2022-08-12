#include "robotoc/mpc/mpc_periodic_configuration_ref.hpp"

#include <stdexcept>
#include <iostream>
#include <cmath>


namespace robotoc {

MPCPeriodicConfigurationRef::MPCPeriodicConfigurationRef(const Eigen::VectorXd& q,
                                                         const double swing_start_time, 
                                                         const double period_active, 
                                                         const double period_inactive,
                                                         const int num_phases_in_period)
  : ConfigurationSpaceRefBase(),
    q_(q), 
    quat_(), 
    has_inactive_contacts_(),
    swing_start_time_(swing_start_time), 
    period_active_(period_active), 
    period_inactive_(period_inactive), 
    period_(period_active+period_inactive),
    num_phases_in_period_(num_phases_in_period) {
  if (period_active < 0.0) {
    throw std::out_of_range(
        "[MPCPeriodicConfigurationRef] invalid argument: period_active must be non-negative!");
  }
  if (period_inactive < 0.0) {
    throw std::out_of_range(
        "[MPCPeriodicConfigurationRef] invalid argument: period_active must be non-negative!");
  }
  if (num_phases_in_period < 1) {
    throw std::out_of_range(
        "[MPCPeriodicConfigurationRef] invalid argument: num_phases_in_period must be positive!");
  }
}


MPCPeriodicConfigurationRef::~MPCPeriodicConfigurationRef() {
}


void MPCPeriodicConfigurationRef::setPeriod(const double swing_start_time, 
                                            const double period_active, 
                                            const double period_inactive,
                                            const int num_phases_in_period) {
  if (period_active < 0.0) {
    throw std::out_of_range(
        "[MPCPeriodicConfigurationRef] invalid argument: period_active must be non-negative!");
  }
  if (period_inactive < 0.0) {
    throw std::out_of_range(
        "[MPCPeriodicConfigurationRef] invalid argument: period_active must be non-negative!");
  }
  if (num_phases_in_period < 1) {
    throw std::out_of_range(
        "[MPCPeriodicConfigurationRef] invalid argument: num_phases_in_period must be positive!");
  }
  swing_start_time_ = swing_start_time;
  period_active_ = period_active;
  period_inactive_ = period_inactive;
  period_ = period_active + period_inactive;
  num_phases_in_period_ = num_phases_in_period;
}


void MPCPeriodicConfigurationRef::setConfigurationRef(
    const std::shared_ptr<ContactSequence>& contact_sequence, 
    const std::shared_ptr<ContactPlannerBase>& foot_step_planner) {
  has_inactive_contacts_.clear();
  for (int phase=0; phase<contact_sequence->numContactPhases(); ++phase) {
    const auto& contact_status = contact_sequence->contactStatus(phase);
    int num_active_contacts = 0;
    for (int i=0; i<contact_status.maxNumContacts(); ++i) {
      if (contact_status.isContactActive(i)) {
        ++num_active_contacts; 
      }
    }
    has_inactive_contacts_.push_back(
        num_active_contacts < contact_status.maxNumContacts());
  }
  const int size = foot_step_planner->size();
  if (quat_.size() < size) {
    quat_.resize(size);
  }
  for (int i=0; i<size; ++i) {
    quat_[i] = Eigen::Quaterniond(foot_step_planner->R(i));
  }
}


void MPCPeriodicConfigurationRef::updateRef(const Robot& robot, 
                                               const GridInfo& grid_info,
                                               Eigen::VectorXd& q_ref) const {
  // if (isActive(grid_info)) { 
  q_ref = q_;
  if ((grid_info.t > swing_start_time_) 
        && has_inactive_contacts_[grid_info.contact_phase]) {
    const int cycle = std::floor((grid_info.t-swing_start_time_)/period_);
    const double rate = (grid_info.t-swing_start_time_-cycle*period_) / period_active_;
    q_ref.template segment<4>(3) 
        = quat_[grid_info.contact_phase].slerp(
              rate, quat_[grid_info.contact_phase+num_phases_in_period_]).coeffs();
  }
  else {
    q_ref.template segment<4>(3) = quat_[grid_info.contact_phase].coeffs();
  }
}


bool MPCPeriodicConfigurationRef::isActive(const GridInfo& grid_info) const {
  // return ((grid_info.t > swing_start_time_) 
  //           && has_inactive_contacts_[grid_info.contact_phase]);
  return true;
}

} // namespace robotoc