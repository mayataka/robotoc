#ifndef ROBOTOC_MPC_FLYING_TROT_HPP_
#define ROBOTOC_MPC_FLYING_TROT_HPP_

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
#include "robotoc/mpc/contact_planner_base.hpp"
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
/// @class MPCFlyingTrot
/// @brief MPC solver for the trot gait of quadrupeds. 
///
class MPCFlyingTrot {
public:
  ///
  /// @brief Construct MPC solver.
  /// @param[in] quadruped_robot Quadruped robot model. 
  /// @param[in] T Length of the horizon. 
  /// @param[in] N Number of the discretization grids of the horizon. 
  /// @param[in] max_steps Maximum number of trot steps over the horizon.
  /// @param[in] nthreads Number of threads used in the parallel computing.
  ///
  MPCFlyingTrot(const Robot& quadruped_robot, const double T, const int N, 
                    const int max_steps, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  MPCFlyingTrot();

  ///
  /// @brief Destructor. 
  ///
  ~MPCFlyingTrot();

  ///
  /// @brief Default copy constructor. 
  ///
  MPCFlyingTrot(const MPCFlyingTrot&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  MPCFlyingTrot& operator=(const MPCFlyingTrot&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MPCFlyingTrot(MPCFlyingTrot&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MPCFlyingTrot& operator=(MPCFlyingTrot&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] foot_step_planner Foot step planner of the gait. 
  /// @param[in] swing_height Swing height of the gait. 
  /// @param[in] flying_time Flying time of the gait. 
  /// @param[in] stance_time Stance time of the gait. 
  /// @param[in] swing_start_time Start time of the gait. 
  ///
  void setGaitPattern(const std::shared_ptr<ContactPlannerBase>& foot_step_planner,
                      const double swing_height, const double flying_time, 
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
  /// @brief Get the solution. 
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
  /// MPCFlyingTrot::updateSolution() must be computed.  
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

  ///
  /// @brief Gets the contact sequence handle.  
  /// @return Shared ptr to the contact sequence.
  ///
  std::shared_ptr<ContactSequence> getContactSequenceHandle();

  ///
  /// @brief Gets the const handle of the MPC solver.  
  /// @return Const reference to the MPC solver.
  ///
  const OCPSolver& getSolver() const { return ocp_solver_; }

  ///
  /// @brief Gets the const handle of the contact sequence.  
  /// @return Const reference to the shared_ptr of the contact sequence.
  ///
  const std::shared_ptr<ContactSequence>& getContactSequence() const { 
    return contact_sequence_; 
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<ContactPlannerBase> foot_step_planner_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  OCPSolver ocp_solver_;
  SolverOptions solver_options_;
  ContactStatus cs_standing_, cs_lfrh_, cs_rflh_, cs_flying_;
  double swing_height_, flying_time_, stance_time_, swing_start_time_, 
         T_, dt_, dtm_, ts_last_, eps_;
  int N_, current_step_, predict_step_;

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

#endif // ROBOTOC_MPC_FLYING_TROT_HPP_