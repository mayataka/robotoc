#include "robotoc/mpc/mpc_periodic_swing_foot_ref.hpp"

#include <stdexcept>
#include <cmath>


namespace robotoc {

MPCPeriodicSwingFootRef::MPCPeriodicSwingFootRef(const int contact_index, 
                                                 const double swing_height, 
                                                 const double swing_start_time, 
                                                 const double period_swing, 
                                                 const double period_stance)
  : TimeVaryingTaskSpace3DRefBase(),
    contact_index_(contact_index),
    is_contact_active_(),
    swing_height_(swing_height),
    swing_start_time_(swing_start_time),
    period_swing_(period_swing), 
    period_stance_(period_stance),
    period_(period_swing+period_stance) {
}


MPCPeriodicSwingFootRef::~MPCPeriodicSwingFootRef() {
}


void MPCPeriodicSwingFootRef::setPeriod(const double swing_start_time, 
                                        const double period_swing, 
                                        const double period_stance) {
  swing_start_time_ = swing_start_time;
  period_swing_ = period_swing;
  period_stance_ = period_stance;
  period_ = period_swing + period_stance;
}


void MPCPeriodicSwingFootRef::setSwingFootRef(
    const std::shared_ptr<ContactSequence>& contact_sequence,
    const std::shared_ptr<FootStepPlannerBase>& foot_step_planner) {
  is_contact_active_.clear();
  for (int phase=0; phase<contact_sequence->numContactPhases(); ++phase) {
    const auto& contact_status = contact_sequence->contactStatus(phase);
    is_contact_active_.push_back(contact_status.isContactActive(contact_index_));
  }
  contact_position_.clear();
  for (const auto& e : foot_step_planner->contactPosition()) {
    contact_position_.push_back(e[contact_index_]);
  }
}


void MPCPeriodicSwingFootRef::update_x3d_ref(const GridInfo& grid_info,
                                             Eigen::VectorXd& x3d_ref) const {
  if (isActive(grid_info)) {
    const int cycle = std::floor((grid_info.t-swing_start_time_)/period_);
    const double rate = (grid_info.t-swing_start_time_-cycle*period_) / period_swing_;
    x3d_ref = (1.0-rate) * contact_position_[grid_info.contact_phase]
                + rate * contact_position_[grid_info.contact_phase+1];
    if (rate < 0.5) {
      x3d_ref.coeffRef(2) += 2 * rate * swing_height_;
    }
    else {
      x3d_ref.coeffRef(2) += 2 * (1-rate) * swing_height_;
    }
  }
}


bool MPCPeriodicSwingFootRef::isActive(const GridInfo& grid_info) const {
  return ((grid_info.t > swing_start_time_) 
            && !is_contact_active_[grid_info.contact_phase]);
}

} // namespace robotoc