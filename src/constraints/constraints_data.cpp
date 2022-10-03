#include "robotoc/constraints/constraints_data.hpp"


namespace robotoc {

ConstraintsData::ConstraintsData(const int time_stage) 
  : is_position_level_valid_(false), 
    is_velocity_level_valid_(false),
    is_acceleration_level_valid_(false),
    is_impact_level_valid_(false) {
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
  else if (time_stage <= -1) {
    is_position_level_valid_     = false;
    is_velocity_level_valid_     = false;
    is_acceleration_level_valid_ = false;
    is_impact_level_valid_      = true;
  }
}


void ConstraintsData::copySlackAndDual(const ConstraintsData& other) {
  if (isPositionLevelValid()) {
    const int size = position_level_data.size();
    for (int i=0; i<size; ++i) {
      position_level_data[i].copySlackAndDual(other.position_level_data[i]);
    }
  }
  if (isVelocityLevelValid()) {
    const int size = velocity_level_data.size();
    for (int i=0; i<size; ++i) {
      velocity_level_data[i].copySlackAndDual(other.velocity_level_data[i]);
    }
  }
  if (isAccelerationLevelValid()) {
    const int size = acceleration_level_data.size();
    for (int i=0; i<size; ++i) {
      acceleration_level_data[i].copySlackAndDual(other.acceleration_level_data[i]);
    }
  }
  if (isImpactLevelValid()) {
    const int size = impact_level_data.size();
    for (int i=0; i<size; ++i) {
      impact_level_data[i].copySlackAndDual(other.impact_level_data[i]);
    }
  }
}

} // namespace robotoc