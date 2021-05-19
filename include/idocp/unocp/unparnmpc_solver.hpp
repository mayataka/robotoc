#ifndef IDOCP_UNPARNMPC_SOLVER_HPP_
#define IDOCP_UNPARNMPC_SOLVER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/utils/aligned_vector.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/unocp/unbackward_correction.hpp"
#include "idocp/unocp/unconstrained_container.hpp"
#include "idocp/line_search/unline_search.hpp"


namespace idocp {

///
/// @class UnParNMPCSolver
/// @brief Optimal control problem solver of unconstrained rigid-body systems 
/// by ParNMPC algorithm. "Unconstrained" means that the system does not have 
/// either a floating base or any contacts.
///
class UnParNMPCSolver {
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
  UnParNMPCSolver(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                  const std::shared_ptr<Constraints>& constraints, 
                  const double T, const int N, const int nthreads=1);

  ///
  /// @brief Default constructor. 
  ///
  UnParNMPCSolver();

  ///
  /// @brief Destructor. 
  ///
  ~UnParNMPCSolver();

  ///
  /// @brief Default copy constructor. 
  ///
  UnParNMPCSolver(const UnParNMPCSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnParNMPCSolver& operator=(const UnParNMPCSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnParNMPCSolver(UnParNMPCSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnParNMPCSolver& operator=(UnParNMPCSolver&&) noexcept = default;

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
  /// @brief Gets the state-feedback gain.
  /// @param[in] stage Time stage of interest. Must be larger than 0 and smaller
  /// than N.
  /// @param[out] Kq The state-feedback gain with respec to the configuration. 
  /// Size must be Robot::dimu() x Robot::dimv().
  /// @param[out] Kv The state-feedback gain with respec to the velocity. 
  /// Size must be Robot::dimu() x Robot::dimv().
  ///
  void getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                            Eigen::MatrixXd& Kv) const;

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
  /// OCPsolver::computeKKTResidual() must be computed.  
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError();

  ///
  /// @return true if the current solution is feasible subject to the 
  /// inequality constraints. Return false if it is not feasible.
  ///
  bool isCurrentSolutionFeasible();

private:
  std::vector<Robot, Eigen::aligned_allocator<Robot>> robots_;
  UnParNMPC parnmpc_;
  UnBackwardCorrection backward_correction_;
  UnLineSearch line_search_;
  UnKKTMatrix unkkt_matrix_;
  UnKKTResidual unkkt_residual_;
  UnSolution s_;
  UnDirection d_;
  int N_, nthreads_;
  double T_, dt_;
  Eigen::VectorXd kkt_error_;

};

} // namespace idocp 


#endif // IDOCP_UNPARNMPC_SOLVER_HPP_ 