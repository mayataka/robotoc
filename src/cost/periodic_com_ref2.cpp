#include "robotoc/cost/periodic_com_ref2.hpp"


namespace robotoc {

PeriodicCoMRef2::PeriodicCoMRef2(const Eigen::Vector3d com_ref0, 
                                 const Eigen::Vector3d vcom_ref, 
                                 const double t0, 
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


PeriodicCoMRef2::~PeriodicCoMRef2() {
}


void PeriodicCoMRef2::update_com_ref(const GridInfo& grid_info,
                                     Eigen::VectorXd& com_ref) const {
  if (is_first_move_half_) {
    if (grid_info.t < t0_) {
      com_ref = com_ref0_;
    }
    else if (grid_info.t < t0_+period_active_) {
      com_ref = com_ref0_ + 0.5 * (grid_info.t-t0_) * vcom_ref_;
    }
    else {
      const int steps = std::floor((grid_info.t-t0_)/period_);
      const double tau = grid_info.t - t0_ - steps * period_;
      if (tau < period_active_) {
        com_ref = com_ref0_ + ((steps-0.5)*period_active_+tau) * vcom_ref_;
      }
      else {
        com_ref = com_ref0_ + (steps+0.5)*period_active_*vcom_ref_;
      }
    }
  }
  else {
    if (grid_info.t < t0_) {
      com_ref = com_ref0_;
    }
    else {
      const int steps = std::floor((grid_info.t-t0_)/period_);
      const double tau = grid_info.t - t0_ - steps * period_;
      if (tau < period_active_) {
        com_ref = com_ref0_ + (steps*period_active_+tau) * vcom_ref_;
      }
      else {
        com_ref = com_ref0_ + (steps+1)*period_active_*vcom_ref_;
      }
    }
  }
}


bool PeriodicCoMRef2::isActive(const GridInfo& grid_info) const {
  return true;
}

} // namespace robotoc