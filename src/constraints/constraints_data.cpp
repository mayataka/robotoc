#include "robotoc/constraints/constraints_data.hpp"


namespace robotoc {

ConstraintsData::ConstraintsData(const int time_stage) 
  : ConstraintsData() {
  setTimeStage(time_stage);
}


ConstraintsData::ConstraintsData()
  : is_position_level_valid_(false), 
    is_velocity_level_valid_(false),
    is_acceleration_level_valid_(false),
    is_impact_level_valid_(false) {
}


void ConstraintsData::setTimeStage(const int time_stage) {
  if (time_stage >= 2) {
    is_position_level_valid_     = true;
    is_velocity_level_valid_     = true;
    is_acceleration_level_valid_ = true;
    is_impact_level_valid_      = false;
  }
  else if (time_stage == 1) {
    is_position_level_valid_     = false;
    is_velocity_level_valid_     = true;
    is_acceleration_level_valid_ = true;
    is_impact_level_valid_      = false;
  }
  else if (time_stage == 0) {
    is_position_level_valid_     = false;
    is_velocity_level_valid_     = false;
    is_acceleration_level_valid_ = true;
    is_impact_level_valid_      = false;
  }
  else {
    is_position_level_valid_     = false;
    is_velocity_level_valid_     = false;
    is_acceleration_level_valid_ = false;
    is_impact_level_valid_      = true;
  }
}

} // namespace robotoc