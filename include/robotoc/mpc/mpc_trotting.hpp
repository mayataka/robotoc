#ifndef ROBOTOC_MPC_TROTTING_HPP_
#define ROBOTOC_MPC_TROTTING_HPP_

#include <memory>

#include "Eigen/Core"
#include "Eigen/Geometry"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/solver/ocp_solver.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/solver/solver_options.hpp"
#include "robotoc/mpc/foot_step_planner_base.hpp"
#include "robotoc/cost/configuration_space_cost.hpp"
#include "robotoc/cost/time_varying_configuration_space_cost.hpp"
#include "robotoc/cost/time_varying_task_space_3d_cost.hpp"
#include "robotoc/cost/time_varying_com_cost.hpp"
#include "robotoc/mpc/mpc_periodic_swing_foot_ref.hpp"
#include "robotoc/mpc/mpc_periodic_com_ref.hpp"
#include "robotoc/mpc/mpc_periodic_configuration_ref.hpp"
#include "robotoc/constraints/joint_position_lower_limit.hpp"
#include "robotoc/constraints/joint_position_upper_limit.hpp"
#include "robotoc/constraints/joint_velocity_lower_limit.hpp"
#include "robotoc/constraints/joint_velocity_upper_limit.hpp"
#include "robotoc/constraints/joint_torques_lower_limit.hpp"
#include "robotoc/constraints/joint_torques_upper_limit.hpp"
#include "robotoc/constraints/friction_cone.hpp"
#include "robotoc/constraints/impulse_friction_cone.hpp"


namespace robotoc {

///
/// @class MPCTrotting
/// @brief MPC solver for the trotting gait of quadrupeds. 
///
class MPCTrotting {
public:
  ///
  /// @brief Construct MPC solver.
  /// @param[in] quadruped_robot Quadruped robot model. 
  /// @param[in] T Length of the horizon. 
  /// @param[in] N Number of the discretization grids of the horizon. 
  /// @param[in] max_steps Maximum number of trotting steps over the horizon.
  /// @param[in] nthreads Number of threads used in the parallel computing.
  ///
  MPCTrotting(const Robot& quadruped_robot, const double T, const int N, 
              const int max_steps, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  MPCTrotting();

  ///
  /// @brief Destructor. 
  ///
  ~MPCTrotting();

  ///
  /// @brief Default copy constructor. 
  ///
  MPCTrotting(const MPCTrotting&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  MPCTrotting& operator=(const MPCTrotting&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MPCTrotting(MPCTrotting&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MPCTrotting& operator=(MPCTrotting&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] foot_step_planner Foot step planner of the gait. 
  /// @param[in] swing_height Swing height of the gait. 
  /// @param[in] swing_time Swing time of the gait. 
  /// @param[in] stance_time Stance time of the gait. 
  /// @param[in] swing_start_time Start time of the gait. 
  ///
  void setGaitPattern(const std::shared_ptr<FootStepPlannerBase>& foot_step_planner,
                      const double swing_height, const double swing_time, 
                      const double stance_time, const double swing_start_time);

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
  /// MPCTrotting::updateSolution() must be computed.  
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
  /// @brief Gets the base rotation cost handle.  
  /// @return Shared ptr to the base rotation cost.
  ///
  std::shared_ptr<TimeVaryingConfigurationSpaceCost> getBaseRotationCostHandle();

  ///
  /// @brief Gets the swing foot task space costs (LF, LH, RF, RH feet) handle.  
  /// @return Shared ptr to the task space cost (LF, LH, RF, RH feet).
  ///
  std::vector<std::shared_ptr<TimeVaryingTaskSpace3DCost>> getSwingFootCostHandle();

  ///
  /// @brief Gets the com cost handle.  
  /// @return Shared ptr to the com cost.
  ///
  std::shared_ptr<TimeVaryingCoMCost> getCoMCostHandle();

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
  std::shared_ptr<FootStepPlannerBase> foot_step_planner_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  OCPSolver ocp_solver_;
  SolverOptions solver_options_;
  ContactStatus cs_standing_, cs_lfrh_, cs_rflh_;
  double swing_height_, swing_time_, stance_time_, swing_start_time_, 
         T_, dt_, dtm_, ts_last_, eps_;
  int N_, current_step_, predict_step_;
  bool enable_stance_phase_;

  std::shared_ptr<ConfigurationSpaceCost> config_cost_;
  std::shared_ptr<TimeVaryingConfigurationSpaceCost> base_rot_cost_;
  std::shared_ptr<TimeVaryingTaskSpace3DCost> LF_foot_cost_, LH_foot_cost_,
                                              RF_foot_cost_, RH_foot_cost_;
  std::shared_ptr<TimeVaryingCoMCost> com_cost_;
  std::shared_ptr<MPCPeriodicConfigurationRef> base_rot_ref_;
  std::shared_ptr<MPCPeriodicSwingFootRef> LF_foot_ref_, LH_foot_ref_,
                                           RF_foot_ref_, RH_foot_ref_;
  std::shared_ptr<MPCPeriodicCoMRef> com_ref_;
  std::shared_ptr<FrictionCone> friction_cone_;

  bool addStep(const double t);

  void resetContactPlacements(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

};

} // namespace robotoc 

#endif // ROBOTOC_MPC_TROTTING_HPP_ 