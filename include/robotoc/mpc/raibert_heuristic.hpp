#ifndef ROBOTOC_RAIBERT_HEURISTIC_HPP_
#define ROBOTOC_RAIBERT_HEURISTIC_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/mpc/contact_planner_base.hpp"


namespace robotoc {

///
/// @class RaibertHeuristic
/// @brief Raibert heuristic for foot step planning. 
///
class RaibertHeuristic {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] stance_time Stance time. 
  /// @param[in] gain Feedback gain of the velocity. 
  ///
  RaibertHeuristic(const double stance_time, const double gain);

  ///
  /// @brief Default constructor. 
  ///
  RaibertHeuristic();

  ///
  /// @brief Destructor. 
  ///
  ~RaibertHeuristic();

  ///
  /// @brief Default copy constructor. 
  ///
  RaibertHeuristic(const RaibertHeuristic&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  RaibertHeuristic& operator=(const RaibertHeuristic&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RaibertHeuristic(RaibertHeuristic&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RaibertHeuristic& operator=(RaibertHeuristic&&) noexcept = default;

  ///
  /// @brief Set parameters.
  /// @param[in] stance_time Stance time. 
  /// @param[in] gain Feedback gain of the velocity. 
  ///
  void setParameters(const double stance_time, const double gain);

  ///
  /// @brief Plans the step length.
  /// @param[in] vcom Current planar velocity of the COM. 
  /// @param[in] vcom_cmd Commanded planar velocity of the COM. 
  /// @param[in] yaw_rate_cmd Commanded yaw-rate of the COM. 
  ///
  void planStepLength(const Eigen::Vector2d& vcom,
                      const Eigen::Vector2d& vcom_cmd, 
                      const double yaw_rate_cmd);

  ///
  /// @brief Gets the step length.
  ///
  const Eigen::Vector3d& stepLength() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  double stance_time_, gain_;
  Eigen::Vector3d step_length_;

};

} // namespace robotoc 

#endif // ROBOTOC_RAIBERT_HEURISTIC_HPP_