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
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  UnconstrParNMPCSolver(const Robot& robot, 
                        const std::shared_ptr<CostFunction>& cost,
                        const std::shared_ptr<Constraints>& constraints, 
                        const double T, const int N, const int nthreads=1);

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
  /// @brief Updates the solution by computing the primal-dual Newon direction.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] line_search If true, filter line search is enabled. If false
  /// filter line search is disabled. Default is false.
  ///
  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, const bool line_search=false);

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
  /// @brief Clear the line search filter. 
  ///
  void clearLineSearchFilter();

  ///
  /// @brief Computes the KKT residual of the optimal control problem. 
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  ///
  void computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v);

  ///
  /// @brief Returns the l2-norm of the KKT residuals.
  /// UnconstrParNMPCsolver::computeKKTResidual() must be computed.  
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError();

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
  Eigen::VectorXd kkt_error_;

};

} // namespace robotoc 

#endif // ROBOTOC_UNCONSTR_PARNMPC_SOLVER_HPP_