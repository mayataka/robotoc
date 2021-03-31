#ifndef IDOCP_CONSTRAINTS_DATA_HPP_
#define IDOCP_CONSTRAINTS_DATA_HPP_

#include <vector>

#include "idocp/constraints/constraint_component_data.hpp"


namespace idocp {

///
/// @class ConstraintsData
/// @brief Data for constraints. Composed of ConstraintComponentData 
/// corrensponding to the components of Constraints.
///
class ConstraintsData {
public:
  ConstraintsData(const int time_stage) {
    if (time_stage >= 2) {
      is_position_level_valid_     = true;
      is_velocity_level_valid_     = true;
      is_acceleration_level_valid_ = true;
      is_impulse_level_valid_      = false;
    }
    else if (time_stage == 1) {
      is_position_level_valid_     = false;
      is_velocity_level_valid_     = true;
      is_acceleration_level_valid_ = true;
      is_impulse_level_valid_      = false;
    }
    else if (time_stage == 0) {
      is_position_level_valid_     = false;
      is_velocity_level_valid_     = false;
      is_acceleration_level_valid_ = true;
      is_impulse_level_valid_      = false;
    }
    else if (time_stage <= -1) {
      is_position_level_valid_     = false;
      is_velocity_level_valid_     = false;
      is_acceleration_level_valid_ = false;
      is_impulse_level_valid_      = true;
    }
  }

  ConstraintsData() 
    : is_position_level_valid_(false), 
      is_velocity_level_valid_(false),
      is_acceleration_level_valid_(false),
      is_impulse_level_valid_(false) {
  }

  ///
  /// @brief Destructor. 
  ///
  ~ConstraintsData() {}

  ///
  /// @brief Default copy constructor. 
  ///
  ConstraintsData(const ConstraintsData&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ConstraintsData& operator=(const ConstraintsData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ConstraintsData(ConstraintsData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ConstraintsData& operator=(ConstraintsData&&) noexcept = default;

  bool isPositionLevelValid() const {
    return is_position_level_valid_;
  }

  bool isVelocityLevelValid() const {
    return is_velocity_level_valid_;
  }

  bool isAccelerationLevelValid() const {
    return is_acceleration_level_valid_;
  }

  bool isImpulseLevelValid() const {
    return is_impulse_level_valid_;
  }

  std::vector<ConstraintComponentData> position_level_data;
  std::vector<ConstraintComponentData> velocity_level_data;
  std::vector<ConstraintComponentData> acceleration_level_data;
  std::vector<ConstraintComponentData> impulse_level_data;

private:
  bool is_position_level_valid_, is_velocity_level_valid_, 
       is_acceleration_level_valid_, is_impulse_level_valid_;

};
  
} // namespace idocp

#endif // IDOCP_CONSTRAINTS_DATA_HPP_