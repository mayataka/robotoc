#ifndef ROBOTOC_GRID_INFO_HPP_
#define ROBOTOC_GRID_INFO_HPP_

#include <iostream>


namespace robotoc {

/// 
/// @enum GridType
/// @brief Type of the grid.
///
enum class GridType {
  Intermediate,
  Impulse,
  Terminal,
};

/// 
/// @class GridInfo
/// @brief Grid information.
///
struct GridInfo {
  ///
  /// @brief Type of the grid.
  ///
  GridType type = GridType::Intermediate;

  ///
  /// @brief Initial time of the horizon.
  ///
  double t0 = 0;

  ///
  /// @brief Time of this grid.
  ///
  double t = 0;

  ///
  /// @brief Time step of this grid.
  ///
  double dt = 0;

  ///
  /// @brief Time step of the next.
  ///
  double dt_next = 0;

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
  /// @brief Flag if the switching constraint is enable or not. 
  ///
  bool switching_constraint = false;

  ///
  /// @brief Sets random. 
  ///
  void setRandom();

  ///
  /// @brief Returns random grid info. 
  ///
  static GridInfo Random();

  ///
  /// @brief Displays the grid info onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const GridInfo& grid_info);

};

} // namespace robotoc

#endif // ROBOTOC_GRID_INFO_HPP_ 