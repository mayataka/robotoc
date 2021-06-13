#ifndef IDOCP_DIRECT_MULTIPLE_SHOOTING_HPP_ 
#define IDOCP_DIRECT_MULTIPLE_SHOOTING_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/utils/aligned_vector.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_time_discretization.hpp"


namespace idocp {

///
/// @class DirectMultipleShooting
/// @brief Direct multiple shooting method of the hybrid optimal control 
/// problems. 
///
class DirectMultipleShooting {
public:
  ///
  /// @brief Construct the direct multiple shooting method.
  /// @param[in] N Number of discretization grids of the horizon. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  DirectMultipleShooting(const int N, const int max_num_impulse, 
                         const int nthreads);

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
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] s Solution. 
  ///
  void initConstraints(OCP& ocp, aligned_vector<Robot>& robots,
                       const ContactSequence& contact_sequence, 
                       const Solution& s) const;

  ///
  /// @brief Computes the KKT residual of optimal control problem in parallel. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTResidual(OCP& ocp, aligned_vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          const Solution& s, KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual) const;

  ///
  /// @brief Computes the KKT system, i.e., the condensed KKT matrix and KKT
  /// residual for Newton's method. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void computeKKTSystem(OCP& ocp, aligned_vector<Robot>& robots,
                        const ContactSequence& contact_sequence,
                        const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                        const Solution& s, KKTMatrix& kkt_matrix, 
                        KKTResidual& kkt_residual) const;

  ///
  /// @brief Returns the l2-norm of the KKT residual of optimal control problem.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_residual KKT residual. 
  ///
  double KKTError(const OCP& ocp, const KKTResidual& kkt_residual);

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
  /// @param[in, out] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integrateSolution(OCP& ocp, const aligned_vector<Robot>& robots,
                         const double primal_step_size,
                         const double dual_step_size,
                         Direction& d, Solution& s) const;

  static const Eigen::VectorXd& q_prev(const OCP& ocp, const Eigen::VectorXd& q, 
                                       const Solution& s, const int time_stage);

private:
  template <typename Algorithm>
  void runParallel(OCP& ocp, aligned_vector<Robot>& robots,
                   const ContactSequence& contact_sequence,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Solution& s, KKTMatrix& kkt_matrix, 
                   KKTResidual& kkt_residual) const;

  int max_num_impulse_, nthreads_;
  Eigen::VectorXd kkt_error_;
};

} // namespace idocp 

#include "idocp/ocp/direct_multiple_shooting.hxx"

#endif // IDOCP_DIRECT_MULTIPLE_SHOOTING_HPP_ 