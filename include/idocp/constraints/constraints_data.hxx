#ifndef IDOCP_CONSTRAINTS_DATA_HXX_
#define IDOCP_CONSTRAINTS_DATA_HXX_

#include "idocp/constraints/constraints_data.hpp"


namespace idocp {

inline ConstraintsData::ConstraintsData(const int time_stage) 
  : is_position_level_valid_(false), 
    is_velocity_level_valid_(false),
    is_acceleration_level_valid_(false),
    is_impulse_level_valid_(false) {
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


inline ConstraintsData::ConstraintsData()
  : is_position_level_valid_(false), 
    is_velocity_level_valid_(false),
    is_acceleration_level_valid_(false),
    is_impulse_level_valid_(false) {
}


inline ConstraintsData::~ConstraintsData() {
}


inline bool ConstraintsData::isPositionLevelValid() const {
  return is_position_level_valid_;
}


inline bool ConstraintsData::isVelocityLevelValid() const {
  return is_velocity_level_valid_;
}


inline bool ConstraintsData::isAccelerationLevelValid() const {
  return is_acceleration_level_valid_;
}


inline bool ConstraintsData::isImpulseLevelValid() const {
  return is_impulse_level_valid_;
}


inline double ConstraintsData::KKTError() const {
  double err = 0.0;
  if (isPositionLevelValid()) {
    for (const auto& data : position_level_data) {
      err += data.KKTError();
    }
  }
  if (isVelocityLevelValid()) {
    for (const auto& data : velocity_level_data) {
      err += data.KKTError();
    }
  }
  if (isAccelerationLevelValid()) {
    for (const auto& data : acceleration_level_data) {
      err += data.KKTError();
    }
  }
  if (isImpulseLevelValid()) {
    for (const auto& data : impulse_level_data) {
      err += data.KKTError();
    }
  }
  return err;
}


inline double ConstraintsData::logBarrier() const {
  double lb = 0.0;
  if (isPositionLevelValid()) {
    for (const auto& data : position_level_data) {
      lb += data.log_barrier;
    }
  }
  if (isVelocityLevelValid()) {
    for (const auto& data : velocity_level_data) {
      lb += data.log_barrier;
    }
  }
  if (isAccelerationLevelValid()) {
    for (const auto& data : acceleration_level_data) {
      lb += data.log_barrier;
    }
  }
  if (isImpulseLevelValid()) {
    for (const auto& data : impulse_level_data) {
      lb += data.log_barrier;
    }
  }
  return lb;
}


inline double ConstraintsData::constraintViolation() const {
  double vio = 0.0;
  if (isPositionLevelValid()) {
    for (const auto& data : position_level_data) {
      vio += data.constraintViolation();
    }
  }
  if (isVelocityLevelValid()) {
    for (const auto& data : velocity_level_data) {
      vio += data.constraintViolation();
    }
  }
  if (isAccelerationLevelValid()) {
    for (const auto& data : acceleration_level_data) {
      vio += data.constraintViolation();
    }
  }
  if (isImpulseLevelValid()) {
    for (const auto& data : impulse_level_data) {
      vio += data.constraintViolation();
    }
  }
  return vio;
}

} // namespace idocp

#endif // IDOCP_CONSTRAINTS_DATA_XXH