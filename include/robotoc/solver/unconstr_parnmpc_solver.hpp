#ifndef ROBOTOC_UNCONSTR_PARNMPC_SOLVER_HPP_
#define ROBOTOC_UNCONSTR_PARNMPC_SOLVER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/unconstr/unconstr_parnmpc.hpp"
#include "robotoc/ocp/solution.hpp"
#include "robotoc/ocp/direction.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/parnmpc/unconstr_backward_correction.hpp"
#include "robotoc/line_search/unconstr_line_search.hpp"
#include "robotoc/solver/solver_options.hpp"
#include "robotoc/solver/solver_statistics.hpp"
#include "robotoc/utils/timer.hpp"


namespace robotoc {

///
/// @class UnconstrParNMPCSolver
/// @brief Optimal control problem solver of unconstrained rigid-body systems 
/// by ParNMPC algorithm. "Unconstrained" means that the system does not have 
/// either a floating base or any contacts.
///
class UnconstrParNMPCSolver {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] parnmpc Optimal control problem. 
  /// @param[in] solver_options Solver options. Default is 
  /// SolverOptions().
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  UnconstrParNMPCSolver(
      const UnconstrParNMPC& parnmpc, 
      const SolverOptions& solver_options=SolverOptions(), 
      const int nthreads=1);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrParNMPCSolver();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrParNMPCSolver();

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrParNMPCSolver(const UnconstrParNMPCSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnconstrParNMPCSolver& operator=(const UnconstrParNMPCSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrParNMPCSolver(UnconstrParNMPCSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrParNMPCSolver& operator=(UnconstrParNMPCSolver&&) noexcept = default;

  ///
  /// @brief Sets the solver option. 
  /// @param[in] solver_options Solver options.  
  ///
  void setSolverOptions(const SolverOptions& solver_options);

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  ///
  void initConstraints();

  ///
  /// @brief Initializes the backward correction solver.
  /// @param[in] t Initial time of the horizon. 
  ///
  void initBackwardCorrection(const double t);

  ///
  /// @brief Performs single Newton-type iteration, computes the primal-dual 
  /// Newon direction, and updates the solution.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  ///
  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  ///
  /// @brief Solves the optimal control problem. Internally calls 
  /// updateSolution().
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] init_solver If true, initializes the solver, that is, calls
  /// initConstraints(), initBackwardCorrection(), and clears the line search 
  /// filter. Default is true.
  ///
  void solve(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v,
             const bool init_solver=true);

  ///
  /// @brief Gets the solver statistics.
  /// @return Solver statistics.
  ///
  const SolverStatistics& getSolverStatistics() const;

  ///
  /// @brief Get the split solution of a time stage. For example, the control 
  /// input torques at the initial stage can be obtained by ocp.getSolution(0).u.
  /// @param[in] stage Time stage of interest. Must be larger than 0 and smaller
  /// than N.
  /// @return Const reference to the split solution of the specified time stage.
  ///
  const SplitSolution& getSolution(const int stage) const;

  ///
  /// @brief Get the solution vector over the horizon. 
  /// @param[in] name Name of the variable. 
  /// @return Solution vector.
  ///
  std::vector<Eigen::VectorXd> getSolution(const std::string& name) const;

  ///
  /// @brief Sets the solution over the horizon. 
  /// @param[in] name Name of the variable. 
  /// @param[in] value Value of the specified variable. 
  ///
  void setSolution(const std::string& name, const Eigen::VectorXd& value);

  ///
  /// @brief Computes the KKT residual of the optimal control problem and 
  /// returns the KKT error, that is, the l2-norm of the KKT residual. 
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @return The KKT error, that is, the l2-norm of the KKT residual.
  ///
  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

  ///
  /// @brief Returns the l2-norm of the KKT residuals using the results of 
  /// UnconstrParNMPCsolver::updateSolution() or UnconstrParNMPCsolver::solve().
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError() const;

  ///
  /// @brief Returns the value of the cost function.
  /// UnconstrParNMPCsolver::updateSolution() or 
  /// UnconstrParNMPCsolver::computeKKTResidual() must be called.  
  /// @return The value of the cost function.
  ///
  double cost() const;

  ///
  /// @return true if the current solution is feasible subject to the 
  /// inequality constraints. Return false if it is not feasible.
  ///
  bool isCurrentSolutionFeasible();

  ///
  ///
  /// @brief Sets a collection of the properties for robot model in this solver. 
  /// @param[in] properties A collection of the properties for the robot model.
  ///
  void setRobotProperties(const RobotProperties& properties);

private:
  aligned_vector<Robot> robots_;
  UnconstrParNMPC parnmpc_;
  UnconstrBackwardCorrection backward_correction_;
  UnconstrLineSearch line_search_;
  KKTMatrix kkt_matrix_;
  KKTResidual kkt_residual_;
  Solution s_;
  Direction d_;
  int N_, nthreads_;
  double T_, dt_;
  SolverOptions solver_options_;
  SolverStatistics solver_statistics_;
  Timer timer_;

};

} // namespace robotoc 

#endif // ROBOTOC_UNCONSTR_PARNMPC_SOLVER_HPP_