#include "robotoc/mpc/raibert_heuristic.hpp"

#include <stdexcept>
#include <iostream>
#include <cassert>
#include <cmath>


namespace robotoc {

RaibertHeuristic::RaibertHeuristic(const double t_stance, const double gain)
  : t_stance_(t_stance),
    gain_(gain),
    step_length_(Eigen::Vector3d::Zero()) {
  try {
    if (t_stance <= 0.0) {
      throw std::out_of_range("invalid argument: t_stance must be positive!");
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


RaibertHeuristic::RaibertHeuristic() {
}


RaibertHeuristic::~RaibertHeuristic() {
}


void RaibertHeuristic::planStepLength(const Eigen::Vector2d& v_com,
                                      const Eigen::Vector2d& v_com_cmd, 
                                      const double yaw_rate_cmd) {
  step_length_.template head<2>() = 0.5 * t_stance_ * v_com
                                    + gain_ * (v_com - v_com_cmd);
}


const Eigen::Vector3d& RaibertHeuristic::stepLength() const {
  return step_length_;
}

} // namespace robotoc 