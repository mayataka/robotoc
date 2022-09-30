#ifndef ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HPP_ 
#define ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HPP_

#include <vector>
#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/core/solution.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/grid_info.hpp"
#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/ocp/intermediate_stage.hpp"
#include "robotoc/ocp/impact_stage.hpp"
#include "robotoc/ocp/terminal_stage.hpp"
#include "robotoc/ocp/ocp_def.hpp"


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
  DirectMultipleShooting(const OCPDef& ocp, const int nthreads);

  ///
  /// @brief Default constructor. 
  ///
  DirectMultipleShooting();

  ///
  /// @brief Default destructor. 
  ///
  ~DirectMultipleShooting() = default;

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
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] s Solution. 
  ///
  void initConstraints(aligned_vector<Robot>& robots,
                       const TimeDiscretization& time_discretization, 
                       const Solution& s);

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] s Solution. 
  ///
  bool isFeasible(aligned_vector<Robot>& robots, 
                  const TimeDiscretization& time_discretization, 
                  const Solution& s);

  ///
  /// @brief Computes the KKT residual of optimal control problem in parallel. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalOCP(aligned_vector<Robot>& robots, 
               const TimeDiscretization& time_discretization, 
               const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
               const Solution& s, KKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual of optimal control problem in parallel. 
  /// @param[in, out] ocp Optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalKKT(aligned_vector<Robot>& robots, 
               const TimeDiscretization& time_discretization, 
               const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
               const Solution& s, KKTMatrix& kkt_matrix, 
               KKTResidual& kkt_residual);

  void computeInitialStateDirection(const Robot& robot,  
                                    const Eigen::VectorXd& q0, 
                                    const Eigen::VectorXd& v0, 
                                    const Solution& s, Direction& d) const;

  PerformanceIndex getEval(const TimeDiscretization& time_discretization) const;

  void computeStepSizes(const TimeDiscretization& time_discretization,
                        Direction& d);

  double maxPrimalStepSize() const;

  double maxDualStepSize() const;

  void integrateSolution(const aligned_vector<Robot>& robots,
                         const TimeDiscretization& time_discretization, 
                         const double primal_step_size,
                         const double dual_step_size,
                         const KKTMatrix& kkt_matrix,
                         Direction& d, Solution& s);

private:
  int nthreads_;
  aligned_vector<OCPData> ocp_data_;
  IntermediateStage intermediate_stage_;
  ImpactStage impact_stage_;
  TerminalStage terminal_stage_;
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;
};

} // namespace robotoc 

#endif // ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HPP_ 