#ifndef IDOCP_UNOCP_SOLVER_HPP_
#define IDOCP_UNOCP_SOLVER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/unocp/unriccati_recursion.hpp"
#include "idocp/unocp/unconstrained_container.hpp"
#include "idocp/line_search/unline_search.hpp"


namespace idocp {

///
/// @class UnOCPSolver
/// @brief Optimal control problem solver by Riccati recursion for 
/// "unconstrained" rigid-body systems. "Unconstrained" means that the system 
/// does not have either a floating base or any contacts.
///
class UnOCPSolver {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  UnOCPSolver(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
              const std::shared_ptr<Constraints>& constraints, const double T, 
              const int N, const int nthreads=1);

  ///
  /// @brief Default constructor. 
  ///
  UnOCPSolver();

  ///
  /// @brief Destructor. 
  ///
  ~UnOCPSolver();

  ///
  /// @brief Default copy constructor. 
  ///
  UnOCPSolver(const UnOCPSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnOCPSolver& operator=(const UnOCPSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnOCPSolver(UnOCPSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnOCPSolver& operator=(UnOCPSolver&&) noexcept = default;

  ///
  /// @brief Initializes the inequality constraints, i.e., set slack variables 
  /// and the Lagrange multipliers of inequality constraints. Based on the 
  /// current solution.
  ///
  void initConstraints();

  ///
  /// @brief Updates the solution by computing the primal-dual Newon direction.
  /// @param[in] t Initial time of the horizon. Current time in MPC. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  /// @param[in] use_line_search If true, filter line search is enabled. If 
  /// false, it is disabled. Default is false.
  ///
  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v, 
                      const bool use_line_search=false);

  ///
  /// @brief Get the const reference to the split solution of a time stage. 
  /// For example, you can get the const reference to the control input torques 
  /// at the initial stage via ocp.getSolution(0).u.
  /// @param[in] stage Time stage of interest. Must be more than 0 and less 
  /// than N.
  /// @return Const reference to the split solution of the specified time stage.
  ///
  const SplitSolution& getSolution(const int stage) const;

  ///
  /// @brief Gets the state-feedback gain for the control input torques.
  /// @param[in] stage Time stage of interest. Must be more than 0 and less 
  /// than N-1.
  /// @param[out] Kq Gain with respec to the configuration. Size must be 
  /// Robot::dimv() x Robot::dimv().
  /// @param[out] Kv Gain with respec to the velocity. Size must be
  /// Robot::dimv() x Robot::dimv().
  ///
  void getStateFeedbackGain(const int stage, Eigen::MatrixXd& Kq, 
                            Eigen::MatrixXd& Kv) const;

  ///
  /// @brief Sets the configuration and velocity over the horizon uniformly. 
  /// @param[in] q Configuration. Size must be Robot::dimq().
  /// @param[in] v Velocity. Size must be Robot::dimv().
  ///
  bool setStateTrajectory(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v);

  ///
  /// @brief Clear the line search filter. 
  ///
  void clearLineSearchFilter();

  ///
  /// @brief Computes the KKT residula of the optimal control problem. 
  /// @param[in] t Current time. 
  /// @param[in] q Initial configuration. Size must be Robot::dimq().
  /// @param[in] v Initial velocity. Size must be Robot::dimv().
  ///
  void computeKKTResidual(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v);

  ///
  /// @brief Returns the squared KKT error norm by using previously computed 
  /// results computed by updateSolution(). The result is not exactly the 
  /// same as the squared KKT error norm of the original optimal control 
  /// problem. The result is the squared norm of the condensed residual. 
  /// However, this variables is sufficiently close to the original KKT error norm.
  /// @return The squared norm of the condensed KKT residual.
  ///
  double KKTError();

  ///
  /// @brief Return true if the current solution is feasible under the 
  /// inequality constraints. Return false if it is not feasible.
  /// @return true if the current solution is feasible under the inequality 
  /// constraints. false if it is not feasible.
  ///
  bool isCurrentSolutionFeasible();

  ///
  /// @brief Creates robot model.
  /// @return Robot model.
  ///
  Robot createRobot() const;

  ///
  /// @brief Get the solution vector. This function is not suitable for 
  /// real-time application, e.g., MPC, since this function reconstructs the 
  /// solution vector object.
  /// @param[in] name Name of the printed variable. 
  ///
  std::vector<Eigen::VectorXd> getSolution(const std::string& name) const;

  ///
  /// @brief Prints the variable into console. 
  /// @param[in] name Name of the printed variable. Default is "all" 
  /// (print all variables).
  /// @param[in] frame_id Index of the end-effector frames. Only used if 
  /// name == "end-effector". Default is {} (do not specify any frames).
  ///
  void printSolution(const std::string& name="all", 
                     const std::vector<int> frames={}) const;

  ///
  /// @brief Save the variable into file. 
  /// @param[in] name Name of the printed variable. 
  /// @param[in] frame_id Index of the end-effector frames. Only used if 
  /// name == "end-effector". Default is {} (do not specify any frames).
  ///
  void saveSolution(const std::string& path_to_file,
                    const std::string& name) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Robot> robots_;
  UnOCP ocp_;
  UnRiccatiRecursion riccati_recursion_;
  UnLineSearch line_search_;
  SplitKKTMatrix terminal_kkt_matrix_;
  SplitKKTResidual terminal_kkt_residual_;
  UnKKTMatrix unkkt_matrix_;
  UnKKTResidual unkkt_residual_;
  UnSolution s_;
  UnDirection d_;
  UnRiccatiFactorization riccati_factorization_;
  int N_, nthreads_;
  double T_, dtau_;
  Eigen::VectorXd primal_step_size_, dual_step_size_, kkt_error_;

};

} // namespace idocp 


#endif // IDOCP_UNOCP_SOLVER_HPP_ 