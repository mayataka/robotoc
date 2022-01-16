#ifndef ROBOTOC_GRID_INFO_HPP_
#define ROBOTOC_GRID_INFO_HPP_

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
};

} // namespace robotoc

#endif // ROBOTOC_GRID_INFO_HPP_ 