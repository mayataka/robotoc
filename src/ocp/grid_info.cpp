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
  phase = i_distr(eng);
  stage = i_distr(eng);
  stage_in_phase = i_distr(eng);
  num_grids_in_phase = i_distr(eng);
}


GridInfo GridInfo::Random() {
  GridInfo grid_info;
  grid_info.setRandom();
  return grid_info;
}


void GridInfo::disp(std::ostream& os) const {
  auto gridTypeToString = [](const GridType& type) {
    switch (type)
    {
    case GridType::Intermediate:
      return "Intermediate";
      break;
    case GridType::Impact:
      return "Impact";
      break;
    case GridType::Lift:
      return "Lift";
      break;
    case GridType::Terminal:
      return "Terminal";
      break;
    default:
      return "";
      break;
    }
  };
  os << "GridInfo: " << "\n";
  os << "  type: " << gridTypeToString(type) << "\n";
  os << "  t0:      " << t0 << "\n";
  os << "  t:       " << t << "\n";
  os << "  dt:      " << dt << "\n";
  os << "  dt_next: " << dt_next << "\n";
  os << "  stage:   " << stage << "\n";
  os << "  phase:   " << phase << "\n";
  os << "  impact_index: " << impact_index << "\n";
  os << "  lift_index:    " << lift_index << "\n";
  os << "  stage_in_phase:     " << stage_in_phase << "\n";
  os << "  num_grids_in_phase: " << num_grids_in_phase << "\n";
  os << "  sto:      " << std::boolalpha << sto << "\n";
  os << "  sto_next: " << std::boolalpha << sto_next << "\n";
  os << "  switching_constraint: " << std::boolalpha << switching_constraint << std::flush;
}


std::ostream& operator<<(std::ostream& os, const GridInfo& grid_info) {
  grid_info.disp(os);
  return os;
}

} // namespace robotoc