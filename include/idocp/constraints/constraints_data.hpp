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

  ///
  /// @brief Checks wheather the position-level constraints are valid or not. 
  /// @return true if the position-level constraints are valid. false otherwise. 
  ///
  bool isPositionLevelValid() const;

  ///
  /// @brief Checks wheather the velocity-level constraints are valid or not. 
  /// @return true if the velocity-level constraints are valid. false otherwise. 
  ///
  bool isVelocityLevelValid() const;

  ///
  /// @brief Checks wheather the acceleration-level constraints are valid or not. 
  /// @return true if the acceleration-level constraints are valid. false 
  /// otherwise. 
  ///
  bool isAccelerationLevelValid() const;

  ///
  /// @brief Checks wheather the impulse-level constraints are valid or not. 
  /// @return true if the impulse-level constraints are valid. false otherwise. 
  ///
  bool isImpulseLevelValid() const;

  ///
  /// @brief Returns the sum of the squared norm of the KKT error 
  /// (primal residual and complementary slackness) of all the constraints. 
  /// @return true if the impulse-level constraints are valid. false otherwise. 
  ///
  double KKTError() const;

  ///
  /// @brief Returns the sum of the log-barrier of the slack variables of all 
  /// the constraints. 
  /// @return The sum of the log-barrier of the slack variables of all 
  /// the constraints. 
  ///
  double logBarrier() const;

  ///
  /// @brief Returns the sum of the l1-norm of the primal violations in the 
  /// constraints.
  /// @return The sum of the l1-norm of the violations in the constraints. 
  ///
  double constraintViolation() const;

  ///
  /// @brief The collection of the position-level constraints data. 
  ///
  std::vector<ConstraintComponentData> position_level_data;

  ///
  /// @brief The collection of the velocity-level constraints data. 
  ///
  std::vector<ConstraintComponentData> velocity_level_data;

  ///
  /// @brief The collection of the acceleration-level constraints data. 
  ///
  std::vector<ConstraintComponentData> acceleration_level_data;

  ///
  /// @brief The collection of the impulse-level constraints data. 
  ///
  std::vector<ConstraintComponentData> impulse_level_data;

private:
  bool is_position_level_valid_, is_velocity_level_valid_, 
       is_acceleration_level_valid_, is_impulse_level_valid_;

};
  
} // namespace idocp

#include "idocp/constraints/constraints_data.hxx"

#endif // IDOCP_CONSTRAINTS_DATA_HPP_