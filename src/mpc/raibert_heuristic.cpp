#include "robotoc/mpc/raibert_heuristic.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>


namespace robotoc {

RaibertHeuristic::RaibertHeuristic(const double period, const double gain)
  : period_(period),
    gain_(gain),
    step_length_(Eigen::Vector3d::Zero()) {
  if (period <= 0.0) {
    throw std::out_of_range("[RaibertHeuristic] invalid argument: period must be positive!");
  }
  if (gain <= 0.0) {
    throw std::out_of_range("[RaibertHeuristic] invalid argument: gain must be positive!");
  }
  if (gain > 1.0) {
    throw std::out_of_range("[RaibertHeuristic] invalid argument: gain must be less than 1.0!");
  }
}


RaibertHeuristic::RaibertHeuristic() 
  : period_(0.0),
    gain_(0.0),
    step_length_(Eigen::Vector3d::Zero()) {
}


RaibertHeuristic::~RaibertHeuristic() {
}


void RaibertHeuristic::setParameters(const double period, const double gain) {
  if (period <= 0.0) {
    throw std::out_of_range("[RaibertHeuristic] invalid argument: period must be positive!");
  }
  if (gain <= 0.0) {
    throw std::out_of_range("[RaibertHeuristic] invalid argument: gain must be positive!");
  }
  period_ = period;
  gain_ = gain;
}


void RaibertHeuristic::planStepLength(const Eigen::Vector2d& vcom,
                                      const Eigen::Vector2d& vcom_cmd, 
                                      const double yaw_rate_cmd) {
  step_length_.template head<2>() = period_ * vcom
                                    + period_ * gain_ * (vcom_cmd - vcom);
}


const Eigen::Vector3d& RaibertHeuristic::stepLength() const {
  return step_length_;
}

} // namespace robotoc 