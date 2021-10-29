#ifndef ROBOTOC_OCP_SOLVER_HPP_
#define ROBOTOC_OCP_SOLVER_HPP_

#include <vector>
#include <memory>
#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/solution.hpp"
#include "robotoc/ocp/direction.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"
#include "robotoc/riccati/riccati_recursion.hpp"
#include "robotoc/line_search/line_search.hpp"


namespace robotoc {

///
/// @class OCPSolver
/// @brief Optimal control problem solver by Riccati recursion. 
///
class OCPSolver {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_sequence Shared ptr to the contact sequence.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCPSolver(const Robot& robot, 
            const std::shared_ptr<ContactSequence>& contact_sequence,
            const std::shared_ptr<CostFunction>& cost,
            const std::shared_ptr<Constraints>& constraints, const double T, 
            const int N, const int nthreads=1);


  ///
  /// @brief Default constructor. 
  ///
  OCPSolver();

  ///
  /// @brief Destructor. 
  ///
  ~OCPSolver();

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
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in] t Initial time of the horizon. 
  ///
  void initConstraints(const double t);

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
  /// @param[in] option Option for the solution. If name == "f" and 
  /// option == "WORLD", the contact forces expressed in the world frame is 
  /// returned. if option is set to other values, these expressed in the local
  /// frame are returned.
  /// @return Solution vector.
  ///
  std::vector<Eigen::VectorXd> getSolution(const std::string& name,
                                           const std::string& option="") const;

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
  /// @brief Extrapolates the solution over the grids on the last contact phase.
  /// Also initializes the slack and dual variables of the inequality 
  /// constraints on such stages.
  /// @param[in] t Initial time of the horizon. 
  ///
  void extrapolateSolutionLastPhase(const double t);

  ///
  /// @brief Extrapolates the solution over the grids on the initial contact 
  /// phase. Also initializes the slack and dual variables of the inequality 
  /// constraints on such stages.
  /// @param[in] t Initial time of the horizon. 
  ///
  void extrapolateSolutionInitialPhase(const double t);

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
  /// OCPsolver::updateSolution() or OCPsolver::computeKKTResidual() must be 
  /// called.  
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError();

  ///
  /// @brief Returns the value of the cost function.
  /// OCPsolver::updateSolution() or OCPsolver::computeKKTResidual() must be 
  /// called.  
  /// @return The value of the cost function.
  ///
  double cost() const;

  ///
  /// @return true if the current solution is feasible subject to the 
  /// inequality constraints. Return false if it is not feasible.
  ///
  bool isCurrentSolutionFeasible();

  ///
  /// @brief Checks wheather the formulation of the discretized optimal control 
  /// problem is tractable or not.
  /// @param[in] t Initial time of the horizon. 
  /// @return true if the optimal control problem is tractable. false if not.
  ///
  bool isFormulationTractable(const double t);

  ///
  /// @brief Checks wheather the switching times are consistent. 
  /// @param[in] t Initial time of the horizon. 
  /// @return true if the switching times are consistent. false if not.
  ///
  bool isSwitchingTimeConsistent(const double t);

  ///
  /// @brief Displays the optimal control problem solver onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const OCPSolver& ocp_solver);

private:
  aligned_vector<Robot> robots_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  DirectMultipleShooting dms_;
  RiccatiRecursion riccati_recursion_;
  LineSearch line_search_;
  OCP ocp_;
  KKTMatrix kkt_matrix_;
  KKTResidual kkt_residual_;
  Solution s_;
  Direction d_;
  RiccatiFactorization riccati_factorization_;

  void discretizeSolution();

};

} // namespace robotoc 

#endif // ROBOTOC_OCP_SOLVER_HPP_ 