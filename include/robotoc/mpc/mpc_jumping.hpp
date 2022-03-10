#ifndef ROBOTOC_MPC_JUMPING_HPP_
#define ROBOTOC_MPC_JUMPING_HPP_

#include <vector>
#include <memory>
#include <limits>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/solver/solver_options.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/robot/se3.hpp"
#include "robotoc/mpc/foot_step_planner_base.hpp"


namespace robotoc {

///
/// @class MPCJumping
/// @brief MPC solver for the jumping control. 
///
class MPCJumping {
public:
  ///
  /// @brief Construct MPC solver.
  /// @param[in] ocp Optimal contro problem. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. 
  ///
  MPCJumping(const OCP& ocp, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  MPCJumping();

  ///
  /// @brief Destructor. 
  ///
  ~MPCJumping();

  ///
  /// @brief Default copy constructor. 
  ///
  MPCJumping(const MPCJumping&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  MPCJumping& operator=(const MPCJumping&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MPCJumping(MPCJumping&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MPCJumping& operator=(MPCJumping&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] foot_step_planner Foot step planner of the jump. 
  /// @param[in] flying_time If STO is enabled, this is initial guess of the 
  /// flying time. Otherwise, this is used as the fixed flying time.
  /// @param[in] min_flying_time Minimum flying time. 
  /// @param[in] ground_time If STO is enabled, this is initial guess of the 
  /// ground time. Otherwise, this is used as the fixed ground time.
  /// @param[in] min_ground_time Minimum time duration after landing. 
  ///
  void setJumpPattern(const std::shared_ptr<FootStepPlannerBase>& foot_step_planner,
                      const double flying_time, const double min_flying_time, 
                      const double ground_time, const double min_ground_time);

  ///
  /// @brief Initializes the optimal control problem solover. 
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] solver_options Solver options for the initialization. 
  /// @param[in] sto If true, the STO algorithm is enabled, that is, the 
  /// lift-off and touch-down timings are optimized. 
  ///
  void init(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
            const SolverOptions& solver_options, const bool sto=false);

  ///
  /// @brief Resets the optimal control problem solver using the previous 
  /// results of init() or reset().
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] solver_options Solver options for the initialization. 
  /// @param[in] sto If true, lift-off and touch-down timings are optimized. 
  ///
  void reset(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             const SolverOptions& solver_options, const bool sto=false);

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
  /// MPCJumping::updateSolution() must be computed.  
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<FootStepPlannerBase> foot_step_planner_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  std::shared_ptr<robotoc::STOConstraints> sto_constraints_;
  OCPSolver ocp_solver_;
  SolverOptions solver_options_;
  ContactStatus cs_ground_, cs_flying_;
  robotoc::Solution s_;
  double flying_time_, min_flying_time_, ground_time_, min_ground_time_,
         T_, dt_, dtm_, t_mpc_start_, eps_;
  int N_, current_step_;

  void resetMinimumDwellTimes(const double t, const double min_dt);

  void resetGoalContactPlacements(const Eigen::VectorXd& q);

  void resetContactPlacements(const Eigen::VectorXd& q);

};

} // namespace robotoc 

#endif // ROBOTOC_MPC_JUMPING_HPP_ 