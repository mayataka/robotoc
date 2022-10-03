#ifndef ROBOTOC_OCP_SOLVER_HPP_
#define ROBOTOC_OCP_SOLVER_HPP_

#include <vector>
#include <memory>
#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/robot_properties.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/core/solution.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"
#include "robotoc/riccati/riccati_recursion.hpp"
#include "robotoc/riccati/riccati_factorization.hpp"
#include "robotoc/line_search/line_search.hpp"
#include "robotoc/line_search/line_search_settings.hpp"
#include "robotoc/sto/switching_time_optimization.hpp"
#include "robotoc/sto/sto_cost_function.hpp"
#include "robotoc/sto/sto_constraints.hpp"
#include "robotoc/solver/solution_interpolator.hpp"
#include "robotoc/solver/solver_options.hpp"
#include "robotoc/solver/solver_statistics.hpp"
#include "robotoc/utils/timer.hpp"


namespace robotoc {

///
/// @class OCPSolver
/// @brief Optimal control problem solver by Riccati recursion. 
///
class OCPSolver {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] ocp Optimal control problem. 
  /// @param[in] solver_options Solver options. Default is SolverOptions().
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  /// @note If you consider the switching time optimization (STO) problem,
  /// please use the other constructor.
  ///
  OCPSolver(const OCP& ocp, 
            const SolverOptions& solver_options=SolverOptions(), 
            const int nthreads=1);

  ///
  /// @brief Default constructor. 
  ///
  OCPSolver();

  ///
  /// @brief Default destructor. 
  ///
  ~OCPSolver() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  OCPSolver(const OCPSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPSolver& operator=(const OCPSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPSolver(OCPSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPSolver& operator=(OCPSolver&&) noexcept = default;

  ///
  /// @brief Sets the solver option. 
  /// @param[in] solver_options Solver options.  
  ///
  void setSolverOptions(const SolverOptions& solver_options);

  ///
  /// @brief Discretizes the problem and reiszes the data structures.
  /// @param[in] t Initial time of the horizon. 
  ///
  void discretize(const double t);

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in] t Initial time of the horizon. 
  ///
  void initConstraints(const double t);

  ///
  /// @brief Performs single Newton-type iteration and updates the solution.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @remark The linear and angular velocities of the floating base are assumed
  /// to be expressed in the body local coordinate.
  ///
  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  ///
  /// @brief Solves the optimal control problem. Internally calls 
  /// updateSolutio() and discretize().
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] init_solver If true, initializes the solver, that is, calls
  /// discretize(), initConstraints(), and clears the line search filter.
  /// Default is true.
  /// @remark The linear and angular velocities of the floating base are assumed
  /// to be expressed in the body local coordinate.
  ///
  void solve(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v,
             const bool init_solver=true);

  ///
  /// @brief Gets the solver statistics.
  /// @return Solver statistics.
  ///
  const SolverStatistics& getSolverStatistics() const;

  ///
  /// @brief Get the solution over the horizon. 
  /// @return const reference to the solution.
  ///
  const Solution& getSolution() const;

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
  /// @param[in] option Option for the solution. If name == "f" and 
  /// option == "WORLD", the contact forces expressed in the world frame is 
  /// returned. if option is set to other values, these expressed in the local
  /// frame are returned.
  /// @return Solution vector.
  ///
  std::vector<Eigen::VectorXd> getSolution(const std::string& name,
                                           const std::string& option="") const;

  ///
  /// @brief Gets of the local LQR policies over the horizon. 
  /// @return const reference to the local LQR policies.
  ///
  const aligned_vector<LQRPolicy>& getLQRPolicy() const;

  ///
  /// @brief Gets the Riccati factorizations. This can be interpreted as 
  /// locally approximated cost-to-go functions. 
  /// @return const reference to the Riccati factorizations.
  ///
  const RiccatiFactorization& getRiccatiFactorization() const;

  ///
  /// @brief Sets the solution guess over the horizon. 
  /// @param[in] s Solution. 
  ///
  void setSolution(const Solution& s);

  ///
  /// @brief Sets the solution guess over the horizon. 
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
  /// @remark The linear and angular velocities of the floating base are assumed
  /// to be expressed in the body local coordinate.
  ///
  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

  ///
  /// @brief Returns the l2-norm of the KKT residuals using the results of 
  /// OCPsolver::updateSolution() or OCPsolver::solve().
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError() const;

  ///
  /// @brief Gets the discretization. 
  /// @return Returns const reference to the time discretization. 
  ///
  const TimeDiscretization& getTimeDiscretization() const;

  ///
  ///
  /// @brief Sets a collection of the properties for robot model in this solver. 
  /// @param[in] properties A collection of the properties for the robot model.
  ///
  void setRobotProperties(const RobotProperties& properties);

  ///
  /// @brief Displays the optimal control problem solver onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const OCPSolver& ocp_solver);

private:
  aligned_vector<Robot> robots_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  TimeDiscretization time_discretization_;
  DirectMultipleShooting dms_;
  SwitchingTimeOptimization sto_;
  RiccatiRecursion riccati_recursion_;
  LineSearch line_search_;
  OCP ocp_;
  KKTMatrix kkt_matrix_;
  KKTResidual kkt_residual_;
  Solution s_;
  Direction d_;
  RiccatiFactorization riccati_factorization_;
  SolutionInterpolator solution_interpolator_;
  SolverOptions solver_options_;
  SolverStatistics solver_statistics_;
  Timer timer_;

  void resizeData();

};

} // namespace robotoc 

#endif // ROBOTOC_OCP_SOLVER_HPP_ 