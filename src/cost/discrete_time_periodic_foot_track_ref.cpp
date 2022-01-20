#include "robotoc/cost/discrete_time_periodic_foot_track_ref.hpp"


namespace robotoc {

DiscreteTimePeriodicFootTrackRef::DiscreteTimePeriodicFootTrackRef(const Eigen::Vector3d x3d0, 
                                             const double step_length, 
                                             const double step_height, 
                                             const int start_phase, 
                                             const int end_phase, 
                                             const int active_phases, 
                                             const int inactive_phases, 
                                             const bool is_first_step_half)
  : TimeVaryingTaskSpace3DRefBase(),
    x3d0_(x3d0),
    step_length_(step_length),
    step_height_(step_height),
    start_phase_(start_phase),
    end_phase_(end_phase),
    active_phases_(active_phases),
    inactive_phases_(inactive_phases),
    is_first_step_half_(is_first_step_half) {
}


DiscreteTimePeriodicFootTrackRef::~DiscreteTimePeriodicFootTrackRef() {
}


void DiscreteTimePeriodicFootTrackRef::update_x3d_ref(const GridInfo& grid_info,
                                           Eigen::VectorXd& x3d_ref) const {
  if (grid_info.contact_phase < start_phase_+active_phases_) {
    const double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                          / static_cast<double>(grid_info.N_phase);
    x3d_ref = x3d0_;
    if (is_first_step_half_) {
      x3d_ref.coeffRef(0) += 0.5 * rate * step_length_;
    }
    else {
      x3d_ref.coeffRef(0) += rate * step_length_;
    }
    if (rate < 0.5) {
      x3d_ref.coeffRef(2) += 2 * rate * step_height_;
    }
    else {
      x3d_ref.coeffRef(2) += 2 * (1-rate) * step_height_;
    }
  }
  else {
    for (int i=1; ; ++i) {
      if (grid_info.contact_phase 
            < start_phase_+i*(active_phases_+inactive_phases_)+active_phases_) {
        x3d_ref = x3d0_;
        const double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                              / static_cast<double>(grid_info.N_phase);
        if (is_first_step_half_) {
          x3d_ref.coeffRef(0) += (i - 0.5 + rate) * step_length_;
        }
        else {
          x3d_ref.coeffRef(0) += (i + rate) * step_length_;
        }
        if (rate < 0.5) {
          x3d_ref.coeffRef(2) += 2 * rate * step_height_;
        }
        else {
          x3d_ref.coeffRef(2) += 2 * (1-rate) * step_height_;
        }
        break;
      }
    }
  }
}


bool DiscreteTimePeriodicFootTrackRef::isActive(const GridInfo& grid_info) const {
  if (grid_info.contact_phase < start_phase_ 
        || grid_info.contact_phase >= end_phase_) {
    return false;
  }
  else {
    for (int cycle=0; ; ++cycle) {
      if (grid_info.contact_phase 
            < start_phase_+cycle*(active_phases_+inactive_phases_)+active_phases_) {
        return true;
      }
      else if (grid_info.contact_phase 
            < start_phase_+(cycle+1)*(active_phases_+inactive_phases_)) {
        return false;
      }
    }
    return false;
  }
}

} // namespace robotoc