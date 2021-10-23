#include "robotoc/cost/periodic_com_ref2.hpp"


namespace robotoc {

PeriodicCoMRef2::PeriodicCoMRef2(const Eigen::Vector3d CoM_ref0, 
                                 const Eigen::Vector3d v_CoM_ref, 
                                 const double t0, 
                                 const double period_active, 
                                 const double period_inactive, 
                                 const bool is_first_move_half)
  : TimeVaryingCoMRefBase(),
    CoM_ref0_(CoM_ref0),
    v_CoM_ref_(v_CoM_ref),
    t0_(t0),
    period_active_(period_active), 
    period_inactive_(period_inactive),
    period_(period_active+period_inactive),
    is_first_move_half_(is_first_move_half) {
}


PeriodicCoMRef2::~PeriodicCoMRef2() {
}


void PeriodicCoMRef2::update_CoM_ref(const double t, 
                                     Eigen::VectorXd& CoM_ref) const {
  if (is_first_move_half_) {
    if (t < t0_) {
      CoM_ref = CoM_ref0_;
    }
    else if (t < t0_+period_active_) {
      CoM_ref = CoM_ref0_ + 0.5 * (t-t0_) * v_CoM_ref_;
    }
    else {
      const int steps = std::floor((t-t0_)/period_);
      const double tau = t - t0_ - steps * period_;
      if (tau < period_active_) {
        CoM_ref = CoM_ref0_ + ((steps-0.5)*period_active_+tau) * v_CoM_ref_;
      }
      else {
        CoM_ref = CoM_ref0_ + (steps+0.5)*period_active_*v_CoM_ref_;
      }
    }
  }
  else {
    if (t < t0_) {
      CoM_ref = CoM_ref0_;
    }
    else {
      const int steps = std::floor((t-t0_)/period_);
      const double tau = t - t0_ - steps * period_;
      if (tau < period_active_) {
        CoM_ref = CoM_ref0_ + (steps*period_active_+tau) * v_CoM_ref_;
      }
      else {
        CoM_ref = CoM_ref0_ + (steps+1)*period_active_*v_CoM_ref_;
      }
    }
  }
}


bool PeriodicCoMRef2::isActive(const double t) const {
  return true;
}

} // namespace robotoc