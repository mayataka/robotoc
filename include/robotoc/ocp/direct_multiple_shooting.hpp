#ifndef ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HPP_ 
#define ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/solution.hpp"
#include "robotoc/ocp/direction.hpp"
#include "robotoc/ocp/kkt_matrix.hpp"
#include "robotoc/ocp/kkt_residual.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/hybrid/time_discretization.hpp"
#include "robotoc/hybrid/grid_info.hpp"


namespace robotoc {

///
/// @class DirectMultipleShooting
/// @brief Direct multiple shooting method of the hybrid optimal control 
/// problems. 
///
class DirectMultipleShooting {
public:
  ///
  /// @brief Construct the direct multiple shooting method.
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  DirectMultipleShooting(const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  DirectMultipleShooting();

  ///
  /// @brief Destructor. 
  ///
  ~DirectMultipleShooting();

  ///
  /// @brief Default copy constructor. 
  ///
  DirectMultipleShooting(const DirectMultipleShooting&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  DirectMultipleShooting& operator=(const DirectMultipleShooting&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  DirectMultipleShooting(DirectMultipleShooting&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  DirectMultipleShooting& operator=(DirectMultipleShooting&&) noexcept = default;

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] s Solution. 
  ///
  bool isFeasible(OCP& ocp, aligned_vector<Robot>& robots,
                  const Solution& s) const;

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] s Solution. 
  ///
  void initConstraints(OCP& ocp, aligned_vector<Robot>& robots,
                       const Solution& s) const;

  ///
  /// @brief Computes the KKT residual of optimal control problem in parallel. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTResidual(OCP& ocp, aligned_vector<Robot>& robots, 
                          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          const Solution& s, KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the KKT system, i.e., the condensed KKT matrix and KKT
  /// residual for Newton's method. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTSystem(OCP& ocp, aligned_vector<Robot>& robots,
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        const Solution& s, KKTMatrix& kkt_matrix, 
                        KKTResidual& kkt_residual) const;

  ///
  /// @brief Returns the l2-norm of the KKT residual of optimal control problem.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_residual KKT residual. 
  ///
  static double KKTError(const OCP& ocp, const KKTResidual& kkt_residual);

  ///
  /// @brief Returns the total value of the cost function.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] include_cost_barrier If true, includes the cost due to the 
  /// barrier function. Default is true.
  ///
  static double totalCost(const OCP& ocp, const bool include_cost_barrier=true);

  ///
  /// @brief Computes the initial state direction.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] q0 Initial configuration. 
  /// @param[in] v0 Initial generalized velocity. 
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  static void computeInitialStateDirection(const OCP& ocp, 
                                           const aligned_vector<Robot>& robots,  
                                           const Eigen::VectorXd& q0, 
                                           const Eigen::VectorXd& v0, 
                                           const Solution& s, Direction& d);

  ///
  /// @brief Integrates the solution in parallel.
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in, out] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integrateSolution(OCP& ocp, const aligned_vector<Robot>& robots,
                         const double primal_step_size,
                         const double dual_step_size,
                         const KKTMatrix& kkt_matrix,
                         Direction& d, Solution& s) const;

  static const Eigen::VectorXd& q_prev(const OCP& ocp, const Eigen::VectorXd& q, 
                                       const Solution& s, const int time_stage);

  static double dts_stage(const OCP& ocp, const Direction& d, const int time_stage);

  static double dts_aux(const OCP& ocp, const Direction& d, const int impulse_index);

  static double dts_lift(const OCP& ocp, const Direction& d, const int lift_index);

private:
  template <typename Algorithm>
  void runParallel(OCP& ocp, aligned_vector<Robot>& robots,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Solution& s, KKTMatrix& kkt_matrix, 
                   KKTResidual& kkt_residual) const;

  int nthreads_;
};

} // namespace robotoc 

#include "robotoc/ocp/direct_multiple_shooting.hxx"

#endif // ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HPP_ 