#include "robotoc/ocp/grid_info.hpp"

#include <random>
#include <cmath>
#include <chrono>


namespace robotoc {

void GridInfo::setRandom() {
  std::random_device rand;
  std::default_random_engine eng(rand());
  std::uniform_real_distribution<double> d_distr(0, 10);
  t0 = d_distr(eng);
  t = t0 + std::abs(d_distr(eng));
  dt = d_distr(eng);
  std::uniform_int_distribution<int> i_distr(0, 100);
  contact_phase = i_distr(eng);
  time_stage = i_distr(eng);
  grid_count_in_phase = i_distr(eng);
  N_phase = i_distr(eng);
}


GridInfo GridInfo::Random() {
  GridInfo grid_info;
  grid_info.setRandom();
  return grid_info;
}


void GridInfo::disp(std::ostream& os) const {
  os << "GridInfo: " << std::endl;
  os << "t0:  " << t0 << std::endl;
  os << "t:  " << t << std::endl;
  os << "dt: " << dt << std::endl;
  os << "contact_phase: " << contact_phase << std::endl;
  os << "time_stage: " << time_stage << std::endl;
  os << "impulse_index: " << impulse_index << std::endl;
  os << "lift_index: " << lift_index << std::endl;
  os << "grid_count_in_phase: " << grid_count_in_phase << std::endl;
  os << "N_phase: " << N_phase << std::endl;
}


std::ostream& operator<<(std::ostream& os, const GridInfo& grid_info) {
  grid_info.disp(os);
  return os;
}

} // namespace robotoc