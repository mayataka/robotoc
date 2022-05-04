#ifndef ROBOTOC_MPC_JUMP_HPP_
#define ROBOTOC_MPC_JUMP_HPP_

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
#include "robotoc/mpc/contact_planner_base.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/constraints/friction_cone.hpp"


namespace robotoc {

///
/// @class MPCJump
/// @brief MPC solver for the jump control. 
///
class MPCJump {
public:
  ///
  /// @brief Construct MPC solver.
  /// @param[in] robot Robot model. 
  /// @param[in] T Length of the horizon. 
  /// @param[in] N Number of the discretization grids of the horizon. 
  /// @param[in] max_steps Maximum number of steps over the horizon.
  /// @param[in] nthreads Number of threads used in the parallel computing.
  ///
  MPCJump(const Robot& robot, const double T, const int N, 
             const int max_steps, const int nthreads);
  ///
  /// @brief Default constructor. 
  ///
  MPCJump();

  ///
  /// @brief Destructor. 
  ///
  ~MPCJump();

  ///
  /// @brief Default copy constructor. 
  ///
  MPCJump(const MPCJump&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  MPCJump& operator=(const MPCJump&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MPCJump(MPCJump&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MPCJump& operator=(MPCJump&&) noexcept = default;

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
  void setJumpPattern(const std::shared_ptr<ContactPlannerBase>& foot_step_planner,
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
  /// @brief Get the solution over the horizon. 
  /// @return const reference to the solution.
  ///
  const Solution& getSolution() const;

  ///
  /// @brief Gets of the local LQR policies over the horizon. 
  /// @return const reference to the local LQR policies.
  ///
  const hybrid_container<LQRPolicy>& getLQRPolicy() const;

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
  /// MPCJump::updateSolution() must be computed.  
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError() const;

  ///
  /// @brief Gets the cost function handle.  
  /// @return Shared ptr to the cost function.
  ///
  std::shared_ptr<CostFunction> getCostHandle();

  ///
  /// @brief Gets the configuration space cost handle.  
  /// @return Shared ptr to the configuration space cost.
  ///
  std::shared_ptr<ConfigurationSpaceCost> getConfigCostHandle();

  ///
  /// @brief Gets the constraints handle.  
  /// @return Shared ptr to the constraints.
  ///
  std::shared_ptr<Constraints> getConstraintsHandle();

  ///
  /// @brief Gets the friction cone constraints handle.  
  /// @return Shared ptr to the friction cone constraints.
  ///
  std::shared_ptr<FrictionCone> getFrictionConeHandle();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<ContactPlannerBase> foot_step_planner_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  std::shared_ptr<robotoc::STOCostFunction> sto_cost_;
  std::shared_ptr<robotoc::STOConstraints> sto_constraints_;
  OCPSolver ocp_solver_;
  SolverOptions solver_options_;
  ContactStatus cs_ground_, cs_flying_;
  robotoc::Solution s_;
  double flying_time_, min_flying_time_, ground_time_, min_ground_time_,
         T_, dt_, dtm_, t_mpc_start_, eps_;
  int N_, current_step_;

  std::shared_ptr<ConfigurationSpaceCost> config_cost_;
  std::shared_ptr<FrictionCone> friction_cone_;

  void resetMinimumDwellTimes(const double t, const double min_dt);

  void resetGoalContactPlacements(const Eigen::VectorXd& q);

  void resetContactPlacements(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

};

} // namespace robotoc 

#endif // ROBOTOC_MPC_JUMP_HPP_