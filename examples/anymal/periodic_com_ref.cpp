#include "periodic_com_ref.hpp"


namespace idocp {

PeriodicCoMRef::PeriodicCoMRef(const Eigen::Vector3d CoM_ref0, 
                               const Eigen::Vector3d v_CoM_ref, const double t0, 
                               const double period_swing, 
                               const double period_stance, bool first_step)
  : TimeVaryingCoMRefBase(),
    CoM_ref0_(CoM_ref0),
    v_CoM_ref_(v_CoM_ref),
    t0_(t0),
    period_swing_(period_swing), 
    period_stance_(period_stance),
    period_(period_swing+period_stance),
    first_step_(first_step) {
}


PeriodicCoMRef::~PeriodicCoMRef() {
}


void PeriodicCoMRef::update_CoM_ref(const double t, 
                                    Eigen::VectorXd& CoM_ref) const {
  if (t < t0_+period_swing_) {
    if (first_step_) 
      CoM_ref = CoM_ref0_ + 0.5 * (t-t0_) * v_CoM_ref_;
    else 
      CoM_ref = CoM_ref0_ + (t-t0_) * v_CoM_ref_;
    return;
  }
  for (int i=1; ; ++i) {
    if (t < t0_+i*period_+period_swing_) {
      const double t1 = (t-t0_-i*period_);
      if (first_step_) 
        CoM_ref = CoM_ref0_ + ((i - 0.5) * period_swing_ + t1) * v_CoM_ref_;
      else 
        CoM_ref = CoM_ref0_ + (i * period_swing_ + t1) * v_CoM_ref_;
      return;
    }
  }
}


bool PeriodicCoMRef::isActive(const double t) const {
  for (int i=0; ; ++i) {
    if (t < t0_+i*period_) {
      return false;
    }
    if (t < t0_+i*period_+period_swing_) {
      return true;
    }
  }
}

} // namespace idocp