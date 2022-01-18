#include "robotoc/cost/periodic_foot_track_ref2.hpp"


namespace robotoc {

PeriodicFootTrackRef2::PeriodicFootTrackRef2(const Eigen::Vector3d x3d0, 
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


PeriodicFootTrackRef2::~PeriodicFootTrackRef2() {
}


void PeriodicFootTrackRef2::update_x3d_ref(const GridInfo& grid_info,
                                           Eigen::VectorXd& x3d_ref) const {
  x3d_ref = x3d0_;
  if (is_first_step_half_) {
    if (grid_info.t < t0_) {
      // do nothing
    }
    else if (grid_info.t < t0_+period_swing_) {
      const double rate = (grid_info.t-t0_) / period_swing_;
      x3d_ref.coeffRef(0) += 0.5 * rate * step_length_;
      if (rate < 0.5) {
        x3d_ref.coeffRef(2) += 2 * rate * step_height_;
      }
      else {
        x3d_ref.coeffRef(2) += 2 * (1-rate) * step_height_;
      }
    }
    else if (grid_info.t < t0_+period_) {
      x3d_ref.coeffRef(0) += 0.5 * step_length_;
    }
    else {
      const int steps = std::floor((grid_info.t-t0_)/period_);
      const double tau = grid_info.t - t0_ - steps * period_;
      if (tau < period_swing_) {
        const double rate = tau / period_swing_;
        x3d_ref.coeffRef(0) += (0.5+(steps-1)+rate) * step_length_;
        if (rate < 0.5) {
          x3d_ref.coeffRef(2) += 2 * rate * step_height_;
        }
        else {
          x3d_ref.coeffRef(2) += 2 * (1-rate) * step_height_;
        }
      }
      else {
        x3d_ref.coeffRef(0) += (0.5+steps) * step_length_;
      }
    }
  }
  else {
    if (grid_info.t < t0_) {
      // do nothing
    }
    else {
      const int steps = std::floor((grid_info.t-t0_)/period_);
      const double tau = grid_info.t - t0_ - steps * period_;
      if (tau < period_swing_) {
        const double rate = tau / period_swing_;
        x3d_ref.coeffRef(0) += (steps+rate) * step_length_;
        if (rate < 0.5) {
          x3d_ref.coeffRef(2) += 2 * rate * step_height_;
        }
        else {
          x3d_ref.coeffRef(2) += 2 * (1-rate) * step_height_;
        }
      }
      else {
        x3d_ref.coeffRef(0) += (steps+1) * step_length_;
      }
    }
  }
}


bool PeriodicFootTrackRef2::isActive(const GridInfo& grid_info) const {
  return true;
}

} // namespace robotoc