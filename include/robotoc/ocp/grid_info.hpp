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
  Impact,
  Lift,
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
  /// @brief Time step of the next grid.
  ///
  double dt_next = 0;

  ///
  /// @brief Phase of this grid.
  ///
  int phase = 0;

  ///
  /// @brief Stage of this grid. 
  ///
  int stage = 0;

  ///
  /// @brief Impulse index of this grid.
  ///
  int impulse_index = -1;

  ///
  /// @brief Lift index of this grid. 
  ///
  int lift_index = -1;

  ///
  /// @brief Stage index counded in the phase. 
  ///
  int stage_in_phase = 0;

  ///
  /// @brief Total number of grids in the phase.
  ///
  int num_grids_in_phase = 0;

  ///
  /// @brief Flag if the STO is enabled in the current phase. 
  ///
  bool sto = false;

  ///
  /// @brief Flag if the STO is enabled in the next phase. 
  ///
  bool sto_next = false;

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