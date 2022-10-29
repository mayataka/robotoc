#ifndef ROBOTOC_CONTROL_POLICY_HPP_
#define ROBOTOC_CONTROL_POLICY_HPP_

#include <iostream>

#include "Eigen/Core"

#include "robotoc/solver/ocp_solver.hpp"


namespace robotoc {

///
/// @class ControlPolicy
/// @brief Control pocily constructed for the MPC solution. 
///
struct ControlPolicy {
  ///
  /// @brief Constructs the policy from the OCP solver. 
  /// @param[in] ocp_solver OCP solver. 
  /// @param[in] t Inquired time of the control. 
  ///
  ControlPolicy(const OCPSolver& ocp_solver, const double t);

  ///
  /// @brief Default constructor. 
  ///
  ControlPolicy();

  ///
  /// @brief Default destructor. 
  ///
  ~ControlPolicy() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  ControlPolicy(const ControlPolicy&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ControlPolicy& operator=(const ControlPolicy&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ControlPolicy(ControlPolicy&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ControlPolicy& operator=(ControlPolicy&&) noexcept = default;

  ///
  /// @brief Sets the control policy from the OCP solver. 
  /// @param[in] ocp_solver OCP solver. 
  /// @param[in] t Inquired time of the control. 
  ///
  void set(const OCPSolver& ocp_solver, const double t);

  ///
  /// @brief Inquired time of the control. 
  ///
  double t;

  ///
  /// @brief Joint torques. Size must be Robot::dimu().
  ///
  Eigen::VectorXd tauJ; 

  ///
  /// @brief Joint positions. Size must be Robot::dimu().
  ///
  Eigen::VectorXd qJ; 

  ///
  /// @brief Joint velocities. Size must be Robot::dimu().
  ///
  Eigen::VectorXd dqJ; 

  ///
  /// @brief Joint position gain. Size must be Robot::dimu() x Robot::dimu().
  ///
  Eigen::MatrixXd Kp; 

  ///
  /// @brief Joint velocity gain. Size must be Robot::dimu() x Robot::dimu().
  ///
  Eigen::MatrixXd Kd; 

  static ControlPolicy FromOCPSolver(const OCPSolver& ocp_solver);

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const ControlPolicy& control_policy);

};

} // namespace robotoc 

#endif // ROBOTOC_CONTROL_POLICY_HPP_