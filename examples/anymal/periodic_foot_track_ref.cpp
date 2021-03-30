#include "periodic_foot_track_ref.hpp"


namespace idocp {

PeriodicFootTrackRef::PeriodicFootTrackRef(const Eigen::Vector3d p0, 
                                           const double step_length, 
                                           const double step_height, 
                                           const double t0, 
                                           const double period_swing, 
                                           const double period_stance, 
                                           bool first_step)
  : TimeVaryingTaskSpace3DRefBase(),
    p0_(p0),
    step_length_(step_length),
    step_height_(step_height),
    t0_(t0),
    period_swing_(period_swing), 
    period_stance_(period_stance),
    period_(period_swing+period_stance),
    first_step_(first_step) {
}


PeriodicFootTrackRef::~PeriodicFootTrackRef() {
}


void PeriodicFootTrackRef::update_q_3d_ref(const double t, 
                                           Eigen::VectorXd& q_3d_ref) const {
  if (t < t0_+period_swing_) {
    const double rate = (t-t0_) / period_swing_;
    q_3d_ref = p0_;
    if (first_step_) 
      q_3d_ref.coeffRef(0) += 0.5 * rate * step_length_;
    else 
      q_3d_ref.coeffRef(0) += rate * step_length_;
    if (rate < 0.5) 
      q_3d_ref.coeffRef(2) += 2 * rate * step_height_;
    else 
      q_3d_ref.coeffRef(2) += 2 * (1-rate) * step_height_;
    return;
  }
  for (int i=1; ; ++i) {
    if (t < t0_+i*period_+period_swing_) {
      q_3d_ref = p0_;
      const double rate = (t-t0_-i*period_) / period_swing_;
      if (first_step_) 
        q_3d_ref.coeffRef(0) += (i - 0.5 + rate) * step_length_;
      else 
        q_3d_ref.coeffRef(0) += (i + rate) * step_length_;
      if (rate < 0.5) 
        q_3d_ref.coeffRef(2) += 2 * rate * step_height_;
      else 
        q_3d_ref.coeffRef(2) += 2 * (1-rate) * step_height_;
      return;
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

} // namespace idocp