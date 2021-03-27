#ifndef IDOCP_OCP_LINEARIZER_HPP_ 
#define IDOCP_OCP_LINEARIZER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/state_constraint_jacobian.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"


namespace idocp {

///
/// @class OCPLinearizer
/// @brief Linearizer of the hybrid optimal control problems. 
///
class OCPLinearizer {
public:
  ///
  /// @brief Construct the linearizer.
  /// @param[in] N Number of discretization grids of the horizon. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  ///
  OCPLinearizer(const int N, const int max_num_impulse, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  OCPLinearizer();

  ///
  /// @brief Destructor. 
  ///
  ~OCPLinearizer();

  ///
  /// @brief Default copy constructor. 
  ///
  OCPLinearizer(const OCPLinearizer&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCPLinearizer& operator=(const OCPLinearizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCPLinearizer(OCPLinearizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCPLinearizer& operator=(OCPLinearizer&&) noexcept = default;

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] s Solution. 
  ///
  void initConstraints(OCP& ocp, std::vector<Robot>& robots,
                       const ContactSequence& contact_sequence, 
                       const Solution& s) const;

  ///
  /// @brief Linearizes the optimal control problem in parallel. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in, out] jac Jacobian of the switching constraints. 
  ///
  void linearizeOCP(OCP& ocp, std::vector<Robot>& robots,
                    const ContactSequence& contact_sequence,
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    const Solution& s, KKTMatrix& kkt_matrix, 
                    KKTResidual& kkt_residual,
                    StateConstraintJacobian& jac) const;

  ///
  /// @brief Computes the KKT residual of optimal control problem in parallel. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in, out] jac Jacobian of the switching constraints. 
  ///
  void computeKKTResidual(OCP& ocp, std::vector<Robot>& robots, 
                          const ContactSequence& contact_sequence,
                          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                          const Solution& s, KKTMatrix& kkt_matrix, 
                          KKTResidual& kkt_residual,
                          StateConstraintJacobian& jac) const;

  ///
  /// @brief Returns the l2-norm of the KKT residual of optimal control problem.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_residual KKT residual. 
  ///
  double KKTError(const OCP& ocp, const KKTResidual& kkt_residual);

  ///
  /// @brief Integrates the solution in parallel.
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots std::vector of Robot.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in, out] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integrateSolution(OCP& ocp, const std::vector<Robot>& robots,
                         const KKTMatrix& kkt_matrix, KKTResidual& kkt_residual,
                         const double primal_step_size,
                         const double dual_step_size,
                         Direction& d, Solution& s) const;

  static const Eigen::VectorXd& q_prev(const OCPDiscretizer& ocp_discretizer, 
                                       const Eigen::VectorXd& q, 
                                       const Solution& s, const int time_stage);


private:
  template <typename Algorithm>
  void runParallel(OCP& ocp, std::vector<Robot>& robots,
                   const ContactSequence& contact_sequence,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Solution& s, KKTMatrix& kkt_matrix, 
                   KKTResidual& kkt_residual,
                   StateConstraintJacobian& jac) const;

  int max_num_impulse_, nthreads_;
  Eigen::VectorXd kkt_error_;
};

} // namespace idocp 

#include "idocp/ocp/ocp_linearizer.hxx"

#endif // IDOCP_OCP_LINEARIZER_HPP_