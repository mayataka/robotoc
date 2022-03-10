#ifndef ROBOTOC_MPC_CRAWLING_HPP_
#define ROBOTOC_MPC_CRAWLING_HPP_

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
#include "robotoc/mpc/crawling_foot_step_planner.hpp"


namespace robotoc {

///
/// @class MPCCrawling
/// @brief MPC solver for the crawling gait of quadrupeds. 
///
class MPCCrawling {
public:
  ///
  /// @brief Construct MPC solver.
  /// @param[in] ocp Optimal contro problem. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. 
  ///
  MPCCrawling(const OCP& ocp, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  MPCCrawling();

  ///
  /// @brief Destructor. 
  ///
  ~MPCCrawling();

  ///
  /// @brief Default copy constructor. 
  ///
  MPCCrawling(const MPCCrawling&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  MPCCrawling& operator=(const MPCCrawling&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  MPCCrawling(MPCCrawling&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  MPCCrawling& operator=(MPCCrawling&&) noexcept = default;

  ///
  /// @brief Sets the gait pattern. 
  /// @param[in] vcom Center-of-mass velocity. 
  /// @param[in] yaw_rate Yaw-rate. 
  /// @param[in] swing_time Swing time of the gait. 
  /// @param[in] initial_lift_time Start time of the gait. 
  ///
  void setGaitPattern(const Eigen::Vector3d& vcom, const double yaw_rate,
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
  /// MPCCrawling::updateSolution() must be computed.  
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError() const;

  ///
  /// @brief Gets the foot step planner handle.
  /// @return Shared ptr to the foot step planner.
  ///
  std::shared_ptr<CrawlingFootStepPlanner> getPlanner();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::shared_ptr<CrawlingFootStepPlanner> foot_step_planner_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  OCPSolver ocp_solver_;
  SolverOptions solver_options_;
  ContactStatus cs_standing_, cs_lf_, cs_lh_, cs_rf_, cs_rh_;
  Eigen::Vector3d vcom_, step_length_;
  double step_height_, swing_time_, initial_lift_time_, 
         t_, T_, dt_, dtm_, ts_last_, eps_;
  int N_, current_step_, predict_step_;

  bool addStep(const double t);

  void resetContactPlacements(const Eigen::VectorXd& q);

};

} // namespace robotoc 

#endif // ROBOTOC_MPC_CRAWLING_HPP_ 