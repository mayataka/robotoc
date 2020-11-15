#ifndef IDOCP_OCP_HPP_
#define IDOCP_OCP_HPP_ 

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
// #include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/constraints/impulse_constraints.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/split_impulse_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorizer.hpp"
#include "idocp/ocp/line_search_filter.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"


namespace idocp {

///
/// @class OCP
/// @brief Optimal control problem solver by Riccati recursion. 
///
class OCP {
public:
  ///
  /// @brief Construct optimal control problem solver.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] num_proc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints, const double T, 
      const int N, const int num_proc=1);

  ///
  /// @brief Default constructor. 
  ///
  OCP();

  ///
  /// @brief Destructor. 
  ///
  ~OCP();

  ///
  /// @brief Default copy constructor. 
  ///
  OCP(const OCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCP& operator=(const OCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCP(OCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCP& operator=(OCP&&) noexcept = default;

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
  bool setStateTrajectory(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

  ///
  /// @brief Sets the configuration and velocity over the horizon by linear 
  //// interpolation. 
  /// @param[in] q0 Initial configuration. Size must be Robot::dimq().
  /// @param[in] v0 Initial velocity. Size must be Robot::dimv().
  /// @param[in] qN Terminal configuration. Size must be Robot::dimq().
  /// @param[in] vN Terminal velocity. Size must be Robot::dimv().
  ///
  bool setStateTrajectory(const Eigen::VectorXd& q0, const Eigen::VectorXd& v0,
                          const Eigen::VectorXd& qN, const Eigen::VectorXd& vN);

  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Sets the contact points over the horizon. 
  /// @param[in] contact_points Contact points over the horizon.
  ///
  void setContactPoint(const std::vector<Eigen::Vector3d>& contact_points);

  ///
  /// @brief Sets the contact points over the horizon by the configuration. 
  /// @param[in] q configuration. Size must be Robot::dimq().
  ///
  void setContactPointByKinematics(const Eigen::VectorXd& q);
  
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
  /// @brief Initializes the inequality constraints, i.e., set slack variables 
  /// and the Lagrange multipliers of inequality constraints. Based on the 
  /// current solution.
  ///
  void initConstraints();

  ///
  /// @brief Prints the solution into console. 
  ///
  void printSolution() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:

  hybrid_container<SplitOCP, SplitImpulseOCP> split_ocps_;
  TerminalOCP terminal_ocp_;
  hybrid_container<SplitSolution, ImpulseSplitSolution> s_;
  hybrid_container<SplitDirection, ImpulseSplitDirection> d_;
  hybrid_container<RiccatiFactorization, RiccatiFactorization> riccati_;
  hybrid_container<StateConstraintRiccatiFactorization, StateConstraintRiccatiFactorization> 
      constraint_factorization_;
  StateConstraintRiccatiFactorizer constraint_factorizer_;
  // hybrid_container<double, double> dtau_;
  std::vector<Robot> robots_;
  ContactSequence contact_sequence_;
  LineSearchFilter filter_;
  double T_, dtau_, step_size_reduction_rate_, min_step_size_;
  int N_, num_proc_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_, costs_, violations_,
                  costs_impulse_, violations_impulse, costs_lift_, 
                  violations_lift_;
  

  void linearizeSplitOCPs(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v);

};

} // namespace idocp 


#endif // IDOCP_OCP_HPP_