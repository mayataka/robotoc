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
  ConstraintsData(const int time_stage);

  ConstraintsData();

  ///
  /// @brief Destructor. 
  ///
  ~ConstraintsData();

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

  bool isPositionLevelValid() const;

  bool isVelocityLevelValid() const;

  bool isAccelerationLevelValid() const;

  bool isImpulseLevelValid() const;

  double KKTError() const;

  double logBarrier() const;

  double constraintViolation() const;

  std::vector<ConstraintComponentData> position_level_data;
  std::vector<ConstraintComponentData> velocity_level_data;
  std::vector<ConstraintComponentData> acceleration_level_data;
  std::vector<ConstraintComponentData> impulse_level_data;

private:
  bool is_position_level_valid_, is_velocity_level_valid_, 
       is_acceleration_level_valid_, is_impulse_level_valid_;

};
  
} // namespace idocp

#include "idocp/constraints/constraints_data.hxx"

#endif // IDOCP_CONSTRAINTS_DATA_HPP_