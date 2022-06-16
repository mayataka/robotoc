#include "robotoc/mpc/raibert_heuristic.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>


namespace robotoc {

RaibertHeuristic::RaibertHeuristic(const double stance_time, const double gain)
  : stance_time_(stance_time),
    gain_(gain),
    step_length_(Eigen::Vector3d::Zero()) {
  try {
    if (stance_time <= 0.0) {
      throw std::out_of_range("invalid argument: stance_time must be positive!");
    }
    if (gain <= 0.0) {
      throw std::out_of_range("invalid argument: gain must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


RaibertHeuristic::RaibertHeuristic() 
  : stance_time_(0.0),
    gain_(0.0),
    step_length_(Eigen::Vector3d::Zero()) {
}


RaibertHeuristic::~RaibertHeuristic() {
}


void RaibertHeuristic::setParameters(const double stance_time, const double gain) {
  try {
    if (stance_time <= 0.0) {
      throw std::out_of_range("invalid argument: stance_time must be positive!");
    }
    if (gain <= 0.0) {
      throw std::out_of_range("invalid argument: gain must be positive!");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  stance_time_ = stance_time;
  gain_ = gain;
}


void RaibertHeuristic::planStepLength(const Eigen::Vector2d& vcom,
                                      const Eigen::Vector2d& vcom_cmd, 
                                      const double yaw_rate_cmd) {
  step_length_.template head<2>() = 0.5 * stance_time_ * vcom
                                    - gain_ * (vcom - vcom_cmd);
}


const Eigen::Vector3d& RaibertHeuristic::stepLength() const {
  return step_length_;
}

} // namespace robotoc 