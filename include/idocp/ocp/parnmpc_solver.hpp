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
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
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
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in] t Initial time of the horizon. 
  ///
  void initConstraints(const double t);

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
  /// @brief Sets the contact status over all of the time stages uniformly. Also, 
  /// disable discrete events over all of the time stages.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatusUniformly(const ContactStatus& contact_status);

  ///
  /// @brief Push back the contact status. Discrete events (impulse and lift)
  /// are also appended to the optimal control problem.
  /// @param[in] contact_status Contact status.
  /// @param[in] switching_time Time of the switch of the contact status.
  ///
  void pushBackContactStatus(const ContactStatus& contact_status, 
                             const double switching_time);

  ///
  /// @brief Sets the contact points to contact statsus with specified contact  
  /// phase. Also set the contact points of the discrete event just before the  
  /// contact phase.
  /// @param[in] contact_phase Contact phase.
  /// @param[in] contact_points Contact points.
  ///
  void setContactPoints(const int contact_phase, 
                        const std::vector<Eigen::Vector3d>& contact_points);

  ///
  /// @brief Pop back a contact status. 
  ///
  void popBackContactStatus();

  ///
  /// @brief Pop front a contact status. 
  ///
  void popFrontContactStatus();

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

  ///
  /// @brief Shows the information of the discretized optimal control problem
  /// onto console.
  ///
  void showInfo() const;

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

  void discretizeSolution();

};

} // namespace idocp 


#endif // IDOCP_PARNMPC_SOLVER_HPP_ 