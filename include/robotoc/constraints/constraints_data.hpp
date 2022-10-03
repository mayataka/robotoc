#ifndef ROBOTOC_CONSTRAINTS_DATA_HPP_
#define ROBOTOC_CONSTRAINTS_DATA_HPP_

#include <vector>

#include "robotoc/constraints/constraint_component_data.hpp"


namespace robotoc {

///
/// @class ConstraintsData
/// @brief Data for constraints. Composed of ConstraintComponentData 
/// corrensponding to the components of Constraints.
///
class ConstraintsData {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] time_stage Time stage. 
  ///
  ///
  ConstraintsData(const int time_stage);

  ///
  /// @brief Default constructor. 
  ///
  ConstraintsData();

  ///
  /// @brief Default destructor. 
  ///
  ~ConstraintsData() = default;

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
/// @brief Sets the time stage info. 
/// @param[in] time_stage Time stage. 
///
  void setTimeStage(const int time_stage);

  ///
  /// @brief Checks wheather the position-level constraints are valid or not. 
  /// @return true if the position-level constraints are valid. false otherwise. 
  ///
  bool isPositionLevelValid() const {
    return is_position_level_valid_;
  }

  ///
  /// @brief Checks wheather the velocity-level constraints are valid or not. 
  /// @return true if the velocity-level constraints are valid. false otherwise. 
  ///
  bool isVelocityLevelValid() const {
    return is_velocity_level_valid_;
  }

  ///
  /// @brief Checks wheather the acceleration-level constraints are valid or not. 
  /// @return true if the acceleration-level constraints are valid. false 
  /// otherwise. 
  ///
  bool isAccelerationLevelValid() const {
    return is_acceleration_level_valid_;
  }

  ///
  /// @brief Checks wheather the impulse-level constraints are valid or not. 
  /// @return true if the impulse-level constraints are valid. false otherwise. 
  ///
  bool isImpulseLevelValid() const {
    return is_impulse_level_valid_;
  }

  ///
  /// @brief Copies the slack and dual variables from another constraint data.
  /// @param[in] other Another constraint data. 
  ///
  void copySlackAndDual(const ConstraintsData& other);

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
  /// @brief Returns the lp norm of the primal feasibility, i.e., the constraint 
  /// violation. Default norm is l1-norm. You can also specify l-infty norm by 
  /// passing Eigen::Infinity as the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the primal feasibility.
  ///
  template <int p=1>
  double primalFeasibility() const;

  ///
  /// @brief Returns the lp norm of the dual feasibility. Default norm is 
  /// l1-norm. You can also specify l-infty norm by passing Eigen::Infinity as 
  /// the template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the dual feasibility.
  ///
  template <int p=1>
  double dualFeasibility() const;

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
  
} // namespace robotoc

#include "robotoc/constraints/constraints_data.hxx"

#endif // ROBOTOC_CONSTRAINTS_DATA_HPP_