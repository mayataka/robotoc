#ifndef ROBOTOC_GRID_INFO_HPP_
#define ROBOTOC_GRID_INFO_HPP_

#include <random>
#include <chrono>
#include <iostream>


namespace robotoc {

/// 
/// @class GridInfo
/// @brief Grid information.
///
struct GridInfo {
  ///
  /// @brief Time of this grid.
  ///
  double t = 0;

  ///
  /// @brief Time step of this grid.
  ///
  double dt = 0;

  ///
  /// @brief Contact phase of this grid.
  ///
  int contact_phase = 0;

  ///
  /// @brief Time stage of this grid. This value is valid and takes 
  /// non-negative value only if this grid is a time stage or termina stage.
  int time_stage = -1;

  ///
  /// @brief Impulse index of this grid. This value is valid and takes 
  /// non-negative value only if this grid is an impulse or aux grid.
  ///
  int impulse_index = -1;

  ///
  /// @brief Lift index of this grid. This value is valid and takes 
  /// non-negative value only if this grid is a lift grid.
  ///
  int lift_index = -1;

  ///
  /// @brief Grid index counded in the contact phase that contains this grid. 
  ///
  int grid_count_in_phase = 0;

  ///
  /// @brief Total number of grids in the contact phase that contains this grid. 
  ///
  int N_phase = 0;

  ///
  /// @brief Sets random. 
  ///
  void setRandom() {
    std::random_device rand;
    std::default_random_engine eng(rand());
    std::uniform_real_distribution<double> d_distr(0, 10);
    t = d_distr(eng);
    dt = d_distr(eng);
    std::uniform_int_distribution<int> i_distr(0, 100);
    contact_phase = i_distr(eng);
    time_stage = i_distr(eng);
    grid_count_in_phase = i_distr(eng);
    N_phase = i_distr(eng);
  }

  ///
  /// @brief Returns random grid info. 
  ///
  static GridInfo Random() {
    GridInfo grid_info;
    grid_info.setRandom();
    return grid_info;
  }

  ///
  /// @brief Displays the grid info onto a ostream.
  ///
  void disp(std::ostream& os) const {
    os << "Grid info: " << std::endl;
    os << "t:  " << t << std::endl;
    os << "dt: " << dt << std::endl;
    os << "contact_phase: " << contact_phase << std::endl;
    os << "time_stage: " << time_stage << std::endl;
    os << "grid_count_in_phase: " << grid_count_in_phase << std::endl;
  }

  friend std::ostream& operator<<(std::ostream& os, const GridInfo& grid_info) {
    grid_info.disp(os);
    return os;
  }

};

} // namespace robotoc

#endif // ROBOTOC_GRID_INFO_HPP_ 