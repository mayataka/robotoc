#include "robotoc/cost/periodic_foot_track_ref.hpp"


namespace robotoc {

PeriodicFootTrackRef::PeriodicFootTrackRef(const Eigen::Vector3d x3d0, 
                                           const double step_length, 
                                           const double step_height, 
                                           const double t0, 
                                           const double period_swing, 
                                           const double period_stance, 
                                           const bool is_first_step_half)
  : TimeVaryingTaskSpace3DRefBase(),
    x3d0_(x3d0),
    step_length_(step_length),
    step_height_(step_height),
    t0_(t0),
    period_swing_(period_swing), 
    period_stance_(period_stance),
    period_(period_swing+period_stance),
    is_first_step_half_(is_first_step_half) {
}


PeriodicFootTrackRef::~PeriodicFootTrackRef() {
}


void PeriodicFootTrackRef::update_x3d_ref(const double t, 
                                          Eigen::VectorXd& x3d_ref) const {
  if (t < t0_+period_swing_) {
    const double rate = (t-t0_) / period_swing_;
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
      if (t < t0_+i*period_+period_swing_) {
        x3d_ref = x3d0_;
        const double rate = (t-t0_-i*period_) / period_swing_;
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


bool PeriodicFootTrackRef::isActive(const double t) const {
  for (int i=0; ; ++i) {
    if (t < t0_+i*period_) {
      return false;
    }
    if (t < t0_+i*period_+period_swing_) {
      return true;
    }
  }
}

} // namespace robotoc