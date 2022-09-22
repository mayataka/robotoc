#ifndef ROBOTOC_SWITCHING_CONSTRAINT_DATA_HPP_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_DATA_HPP_

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/switching_constraint_residual.hpp"
#include "robotoc/core/switching_constraint_jacobian.hpp"


namespace robotoc {

///
/// @class SwitchingConstraintData
/// @brief Data for the switching constraint.
///
class SwitchingConstraintData {
public:
  ///
  /// @brief Constructs a switching constraint.
  /// @param[in] robot Robot model. 
  ///
  SwitchingConstraintData(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  SwitchingConstraintData();

  ///
  /// @brief Default destructor. 
  ///
  ~SwitchingConstraintData() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingConstraintData(const SwitchingConstraintData&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SwitchingConstraintData& operator=(const SwitchingConstraintData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingConstraintData(SwitchingConstraintData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingConstraintData& operator=(SwitchingConstraintData&&) noexcept = default;

  Eigen::VectorXd q;

  Eigen::VectorXd dq;

  Eigen::VectorXd PqT_xi;
};

} // namespace robotoc

#endif // ROBOTOC_SWITCHING_CONSTRAINT_DATA_HPP_