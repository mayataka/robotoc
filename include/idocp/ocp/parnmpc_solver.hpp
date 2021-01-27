#ifndef IDOCP_PARNMPC_SOLVER_HPP_
#define IDOCP_PARNMPC_SOLVER_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/parnmpc_linearizer.hpp"
#include "idocp/ocp/backward_correction_solver.hpp"
#include "idocp/line_search/line_search.hpp"


namespace idocp {

///
/// @class ParNMPCSolver
/// @brief Optimal control problem solver by ParNMPC algorithm. 
///
class ParNMPCSolver {
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
  ParNMPCSolver(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
                const std::shared_ptr<Constraints>& constraints, const double T, 
                const int N, const int max_num_impulse=0, const int nthreads=1);

  ///
  /// @brief Default constructor. 
  ///
  ParNMPCSolver();

  ///
  /// @brief Destructor. 
  ///
  ~ParNMPCSolver();

  ///
  /// @brief Default copy constructor. 
  ///
  ParNMPCSolver(const ParNMPCSolver&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ParNMPCSolver& operator=(const ParNMPCSolver&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ParNMPCSolver(ParNMPCSolver&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ParNMPCSolver& operator=(ParNMPCSolver&&) noexcept = default;

  ///
  /// @brief Initializes the inequality constraints, i.e., set slack variables 
  /// and the Lagrange multipliers of inequality constraints. Based on the 
  /// current solution.
  ///
  void initConstraints();

  void initBackwardCorrection(const double t);

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

  void setSolution(const std::string& name, const Eigen::VectorXd& value);

  void setContactStatusUniformly(const ContactStatus& contact_status);

  void pushBackContactStatus(const ContactStatus& contact_status, 
                             const double switching_time,
                             const double t);

  void setContactPoints(const int contact_phase, 
                        const std::vector<Eigen::Vector3d>& contact_points);

  ///
  /// @brief Pop back the discrete event. Contact status after discrete event 
  /// is also removed. 
  ///
  void popBackDiscreteEvent();

  ///
  /// @brief Pop front the discrete event. Contact status before the front 
  /// discrete event is also removed. 
  ///
  void popFrontDiscreteEvent();

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
  ContactSequence contact_sequence_;
  ParNMPCLinearizer parnmpc_linearizer_;
  BackwardCorrectionSolver backward_correction_solver_;
  LineSearch line_search_;
  ParNMPC parnmpc_;
  BackwardCorrection backward_correction_;
  KKTMatrix kkt_matrix_;
  KKTResidual kkt_residual_;
  Solution s_;
  Direction d_;
  int N_, nthreads_;

  void discretizeSolution();

};

} // namespace idocp 


#endif // IDOCP_PARNMPC_SOLVER_HPP_ 