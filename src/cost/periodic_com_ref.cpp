#include "idocp/cost/periodic_com_ref.hpp"


namespace idocp {

PeriodicCoMRef::PeriodicCoMRef(const Eigen::Vector3d CoM_ref0, 
                               const Eigen::Vector3d v_CoM_ref, const double t0, 
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


PeriodicCoMRef::~PeriodicCoMRef() {
}


void PeriodicCoMRef::update_CoM_ref(const double t, 
                                    Eigen::VectorXd& CoM_ref) const {
  if (t < t0_+period_active_) {
    if (is_first_move_half_) {
      CoM_ref = CoM_ref0_ + 0.5 * (t-t0_) * v_CoM_ref_;
    }
    else {
      CoM_ref = CoM_ref0_ + (t-t0_) * v_CoM_ref_;
    }
  }
  else {
    for (int i=1; ; ++i) {
      if (t < t0_+i*period_+period_active_) {
        const double t1 = (t-t0_-i*period_);
        if (is_first_move_half_) {
          CoM_ref = CoM_ref0_ + ((i-0.5)*period_active_+t1) * v_CoM_ref_;
        } 
        else {
          CoM_ref = CoM_ref0_ + (i*period_active_+t1) * v_CoM_ref_;
        }
        break;
      }
    }
  }
}


bool PeriodicCoMRef::isActive(const double t) const {
  for (int i=0; ; ++i) {
    if (t < t0_+i*period_) {
      return false;
    }
    if (t < t0_+i*period_+period_active_) {
      return true;
    }
  }
}

} // namespace idocp