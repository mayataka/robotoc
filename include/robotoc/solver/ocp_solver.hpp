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
#include "robotoc/line_search/line_search_settings.hpp"
#include "robotoc/hybrid/switching_time_optimization.hpp"
#include "robotoc/hybrid/sto_cost_function.hpp"
#include "robotoc/hybrid/sto_constraints.hpp"
#include "robotoc/solver/solver_options.hpp"


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
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] solver_options Solver options. Default is SolverOptions::defaultOptions().
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  /// @note If you consider the switching time optimization (STO) problem,
  /// please use the other constructor.
  ///
  OCPSolver(const OCP& ocp, 
            const std::shared_ptr<ContactSequence>& contact_sequence, 
            const SolverOptions& solver_options=SolverOptions::defaultOptions(), 
            const int nthreads=1);

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
  /// @brief Sets the solver option. 
  /// @param[in] solver_options Solver options.  
  ///
  void setSolverOptions(const SolverOptions& solver_options);

  ///
  /// @brief Applies mesh refinement if the discretization method is   
  /// DiscretizationMethod::PhaseBased. Also initializes the constraints 
  /// if the mesh refiement is carried out.
  /// @param[in] t Initial time of the horizon. 
  ///
  void meshRefinement(const double t);

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
  ///
  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  ///
  /// @brief Solves the optimal control problem. Internally calls 
  /// updateSolutio() and meshRefinement().
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] init_solver If true, initializes the solver, that is, calls
  /// meshRefinement(), initConstraints(), and clears the line search filter.
  /// Default is true.
  ///
  void solve(const double t, const Eigen::VectorXd& q, const Eigen::VectorXd& v,
             const bool init_solver=true);

  ///
  /// @brief Get the split solution of a time stage. For example, the control 
  /// input torques at the initial stage can be obtained by ocp.getSolution(0).u.
  /// @param[in] stage Time stage of interest. Must be larger than 0 and smaller
  /// than N.
  /// @return Const reference to the split solution of the specified time stage.
  ///
  const SplitSolution& getSolution(const int stage) const;

  ///
  /// @brief Get the solution over the horizon. 
  /// @return const reference to the solution.
  ///
  const Solution& getSolution() const;

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
  /// OCPsolver::updateSolution() or OCPsolver::solve().
  /// @return The l2-norm of the KKT residual.
  ///
  double KKTError() const;

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
  bool isCurrentSolutionFeasible(const bool verbose=false);

  ///
  /// @brief OCP discretization. 
  /// @return Returns const reference to the internal OCP discretization. 
  ///
  const HybridOCPDiscretization& getOCPDiscretization() const;

  ///
  /// @brief Set settings for line search. Defalt is 
  /// LineSearchSettings::defaultSettings().
  /// @param[in] settings Line search settings.
  ///
  void setLineSearchSettings(
      const LineSearchSettings& settings=LineSearchSettings::defaultSettings());

  ///
  /// @brief Displays the optimal control problem solver onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const OCPSolver& ocp_solver);

private:
  aligned_vector<Robot> robots_;
  std::shared_ptr<ContactSequence> contact_sequence_;
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
  SolverOptions solver_options_;

  void discretizeSolution();

};

} // namespace robotoc 

#endif // ROBOTOC_OCP_SOLVER_HPP_ 