#ifndef ROBOTOC_MPC_QUADRUPEDAL_TROTTING_HPP_
#define ROBOTOC_MPC_QUADRUPEDAL_TROTTING_HPP_

#include <vector>
#include <memory>
#include <limits>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/solver/solver_options.hpp"


namespace robotoc {

///
/// @class MPCQuadrupedalTrotting
/// @brief MPC solver for the trotting gait of quadrupeds. 
///
class MPCQuadrupedalTrotting {
public:
  ///
  /// @brief Construct MPC solver.
  /// @param[in] ocp Optimal contro problem. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. 
  ///
  MPCQuadrupedalTrotting(const OCP& ocp, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  MPCQuadrupedalTrotting();

  ///
  /// @brief Destructor. 
  ///
  ~MPCQuadrupedalTrotting();

  ///
  /// @brief Default copy constructor. 
  ///
  MPCQuadrupedalTrotting(const MPCQuadrupedalTrotting&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  MPCQuadrupedalTrotting& operator=(const MPCQuadrupedalTrotting&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MPCQuadrupedalTrotting(MPCQuadrupedalTrotting&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MPCQuadrupedalTrotting& operator=(MPCQuadrupedalTrotting&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] vcom_cmd Center-of-mass velocity command. 
  /// @param[in] yaw_rate_cmd Yaw-rate command. 
  /// @param[in] swing_time Swing time of the gait. 
  /// @param[in] initial_lift_time Start time of the gait. 
  ///
  void setGaitPattern(const Eigen::Vector3d& vcom_cmd, const double yaw_rate_cmd,
                      const double swing_time, const double initial_lift_time);

  ///
  /// @brief Initializes the optimal control problem solover. 
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] solver_options Solver options for the initialization. 
  ///
  void init(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const SolverOptions& solver_options);

  ///
  /// @brief Sets the solver options. 
  /// @param[in] solver_options Solver options.  
  ///
  void setSolverOptions(const SolverOptions& solver_options);

  ///
  /// @brief Updates the solution by iterationg the Newton-type method.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] dt Sampling time of MPC. Must be positive.
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Velocity. Size must be Robot::dimv().
  ///
  void updateSolution(const double t, const double dt, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  ///
  /// @brief Get the initial control input.
  /// @return Const reference to the control input.
  ///
  const Eigen::VectorXd& getInitialControlInput() const;

  ///
  /// @brief Computes the KKT residual of the optimal control problem. 
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  ///
  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

  ///
  /// @brief Returns the l2-norm of the KKT residuals.
  /// MPCQuadrupedalTrotting::updateSolution() must be computed.  
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Robot robot_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  OCPSolver ocp_solver_;
  SolverOptions solver_options_;
  ContactStatus cs_standing_, cs_lfrh_, cs_rflh_;
  std::vector<Eigen::Vector3d> contact_positions_, contact_positions_curr_, 
                               contact_positions_prev_;
  Eigen::Vector3d vcom_cmd_, step_length_, com_, com_curr_, com_prev_;
  Eigen::Matrix3d R_, R_yaw_cmd_;
  double step_height_, swing_time_, initial_lift_time_, 
         t_, T_, dt_, dtm_, ts_last_, eps_;
  int N_, current_step_, predict_step_;

  bool addStep(const double t);

  void resetContactPlacements(const Eigen::VectorXd& q);

};

} // namespace robotoc 

#endif // ROBOTOC_MPC_QUADRUPEDAL_TROTTING_HPP_ 