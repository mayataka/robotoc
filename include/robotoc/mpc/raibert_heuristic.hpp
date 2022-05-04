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
#include "robotoc/mpc/foot_step_planner_base.hpp"


namespace robotoc {

///
/// @class RaibertHeuristic
/// @brief Raibert heuristic for foot step planning. 
///
class RaibertHeuristic {
public:
  ///
  /// @brief Constructs the planner.
  /// @param[in] t_stance Stance time. 
  /// @param[in] gain Feedback gain of the velocity. 
  ///
  RaibertHeuristic(const double t_stance, const double gain);

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
  /// @param[in] t_stance Stance time. 
  /// @param[in] gain Feedback gain of the velocity. 
  ///
  void setParameters(const double t_stance, const double gain);

  ///
  /// @brief Plans the step length.
  /// @param[in] v_com Current planar velocity of the COM. 
  /// @param[in] v_com_cmd Commanded planar velocity of the COM. 
  /// @param[in] yaw_rate_cmd Commanded yaw-rate of the COM. 
  ///
  void planStepLength(const Eigen::Vector2d& v_com,
                      const Eigen::Vector2d& v_com_cmd, 
                      const double yaw_rate_cmd);

  ///
  /// @brief Gets the step length.
  ///
  const Eigen::Vector3d& stepLength() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  double t_stance_, gain_;
  Eigen::Vector3d step_length_;

};

} // namespace robotoc 

#endif // ROBOTOC_RAIBERT_HEURISTIC_HPP_