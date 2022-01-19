#include "robotoc/cost/periodic_com_ref2.hpp"


namespace robotoc {

PeriodicCoMRef2::PeriodicCoMRef2(const Eigen::Vector3d com_ref0, 
                                 const Eigen::Vector3d com_step, 
                                 const int start_phase, const int end_phase, 
                                 const int active_phases, 
                                 const int inactive_phases, 
                                 const bool is_first_move_half)
  : TimeVaryingCoMRefBase(),
    com_ref0_(com_ref0),
    com_step_(com_step),
    start_phase_(start_phase),
    end_phase_(end_phase),
    active_phases_(active_phases),
    inactive_phases_(inactive_phases),
    is_first_move_half_(is_first_move_half) {
}


PeriodicCoMRef2::~PeriodicCoMRef2() {
}


void PeriodicCoMRef2::update_com_ref(const GridInfo& grid_info,
                                     Eigen::VectorXd& com_ref) const {
  if (grid_info.contact_phase < start_phase_+active_phases_) {
    const double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                          / static_cast<double>(grid_info.N_phase);
    if (is_first_move_half_) {
      com_ref = com_ref0_ + 0.5 * rate * com_step_;
    }
    else {
      com_ref = com_ref0_ + rate * com_step_;
    }
  }
  else {
    for (int i=1; ; ++i) {
      if (grid_info.contact_phase 
            < start_phase_+i*(active_phases_+inactive_phases_)+active_phases_) {
        const double rate = static_cast<double>(grid_info.grid_count_in_phase) 
                              / static_cast<double>(grid_info.N_phase);
        if (is_first_move_half_) {
          com_ref = com_ref0_ + ((i-0.5)+rate) * com_step_;
        } 
        else {
          com_ref = com_ref0_ + (i+rate) * com_step_;
        }
        break;
      }
    }
  }
}


bool PeriodicCoMRef2::isActive(const GridInfo& grid_info) const {
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