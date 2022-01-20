#include "robotoc/cost/periodic_com_ref.hpp"


namespace robotoc {

PeriodicCoMRef::PeriodicCoMRef(const Eigen::Vector3d com_ref0, 
                               const Eigen::Vector3d vcom_ref, const double t0, 
                               const double period_active, 
                               const double period_inactive, 
                               const bool is_first_move_half)
  : TimeVaryingCoMRefBase(),
    com_ref0_(com_ref0),
    vcom_ref_(vcom_ref),
    t0_(t0),
    period_active_(period_active), 
    period_inactive_(period_inactive),
    period_(period_active+period_inactive),
    is_first_move_half_(is_first_move_half) {
}


PeriodicCoMRef::~PeriodicCoMRef() {
}


void PeriodicCoMRef::update_com_ref(const GridInfo& grid_info,
                                    Eigen::VectorXd& com_ref) const {
  if (grid_info.t < t0_+period_active_) {
    if (is_first_move_half_) {
      com_ref = com_ref0_ + 0.5 * (grid_info.t-t0_) * vcom_ref_;
    }
    else {
      com_ref = com_ref0_ + (grid_info.t-t0_) * vcom_ref_;
    }
  }
  else {
    for (int i=1; ; ++i) {
      if (grid_info.t < t0_+i*period_+period_active_) {
        const double t1 = (grid_info.t-t0_-i*period_);
        if (is_first_move_half_) {
          com_ref = com_ref0_ + ((i-0.5)*period_active_+t1) * vcom_ref_;
        } 
        else {
          com_ref = com_ref0_ + (i*period_active_+t1) * vcom_ref_;
        }
        break;
      }
    }
  }
}


bool PeriodicCoMRef::isActive(const GridInfo& grid_info) const {
  for (int i=0; ; ++i) {
    if (grid_info.t < t0_+i*period_) {
      return false;
    }
    if (grid_info.t < t0_+i*period_+period_active_) {
      return true;
    }
  }
}

} // namespace robotoc