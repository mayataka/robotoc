#ifndef ROBOTOC_CONSTRAINTS_DATA_HXX_
#define ROBOTOC_CONSTRAINTS_DATA_HXX_

#include "robotoc/constraints/constraints_data.hpp"


namespace robotoc {

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
  if (isImpactLevelValid()) {
    for (const auto& data : impact_level_data) {
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
  if (isImpactLevelValid()) {
    for (const auto& data : impact_level_data) {
      lb += data.log_barrier;
    }
  }
  return lb;
}


template <int p>
inline double ConstraintsData::primalFeasibility() const {
  double feasibility = 0.0;
  if (isPositionLevelValid()) {
    for (const auto& data : position_level_data) {
      feasibility += data.template primalFeasibility<p>();
    }
  }
  if (isVelocityLevelValid()) {
    for (const auto& data : velocity_level_data) {
      feasibility += data.template primalFeasibility<p>();
    }
  }
  if (isAccelerationLevelValid()) {
    for (const auto& data : acceleration_level_data) {
      feasibility += data.template primalFeasibility<p>();
    }
  }
  if (isImpactLevelValid()) {
    for (const auto& data : impact_level_data) {
      feasibility += data.template primalFeasibility<p>();
    }
  }
  return feasibility;
}


template <int p>
inline double ConstraintsData::dualFeasibility() const {
  double feasibility = 0.0;
  if (isPositionLevelValid()) {
    for (const auto& data : position_level_data) {
      feasibility += data.template dualFeasibility<p>();
    }
  }
  if (isVelocityLevelValid()) {
    for (const auto& data : velocity_level_data) {
      feasibility += data.template dualFeasibility<p>();
    }
  }
  if (isAccelerationLevelValid()) {
    for (const auto& data : acceleration_level_data) {
      feasibility += data.template dualFeasibility<p>();
    }
  }
  if (isImpactLevelValid()) {
    for (const auto& data : impact_level_data) {
      feasibility += data.template dualFeasibility<p>();
    }
  }
  return feasibility;
}

} // namespace robotoc

#endif // ROBOTOC_CONSTRAINTS_DATA_XXH