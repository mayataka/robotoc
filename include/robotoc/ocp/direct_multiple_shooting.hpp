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
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/grid_info.hpp"
#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/ocp/intermediate_stage.hpp"
#include "robotoc/ocp/impact_stage.hpp"
#include "robotoc/ocp/terminal_stage.hpp"


namespace robotoc {

///
/// @class DirectMultipleShooting
/// @brief Direct multiple shooting method of the optimal control 
/// problems. 
///
class DirectMultipleShooting {
public:
  ///
  /// @brief Construct the direct multiple shooting method.
  /// @param[in] ocp Optimal control problem. 
  /// @param[in] nthreads Number of the threads of the parallel computations.
  /// Must be positive. 
  ///
  DirectMultipleShooting(const OCP& ocp, const int nthreads);

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
  /// @brief Sets the number of threads of the parallel computations.
  /// @param[in] nthreads Number of the threads of the parallel computations.
  /// Must be positive. 
  ///
  void setNumThreads(const int nthreads);

  ///
  /// @brief Initializes the priaml-dual interior point method for inequality 
  /// constraints. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] s Solution. 
  ///
  void initConstraints(aligned_vector<Robot>& robots,
                       const TimeDiscretization& time_discretization, 
                       const Solution& s);

  ///
  /// @brief Checks whether the solution is feasible under inequality constraints.
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] s Solution. 
  ///
  bool isFeasible(aligned_vector<Robot>& robots, 
                  const TimeDiscretization& time_discretization, 
                  const Solution& s);

  ///
  /// @brief Computes the cost and constraint violations. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] kkt_residual KKT residual. 
  ///
  void evalOCP(aligned_vector<Robot>& robots, 
               const TimeDiscretization& time_discretization, 
               const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
               const Solution& s, KKTResidual& kkt_residual);

  ///
  /// @brief Computes the KKT residual and matrix. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
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

  ///
  /// @brief Computes the initial state direction. 
  /// @param[in] robot Robot model.
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in, out] d Direction. 
  ///
  void computeInitialStateDirection(const Robot& robot,  
                                    const Eigen::VectorXd& q, 
                                    const Eigen::VectorXd& v, 
                                    const Solution& s, Direction& d) const;

  ///
  /// @brief Gets the performance index of the evaluation. 
  /// @return const reference to the performance index.
  ///
  const PerformanceIndex& getEval() const;

  ///
  /// @brief Computes the step sizes via the fraction-to-boundary-rule.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in, out] d Direction. 
  ///
  void computeStepSizes(const TimeDiscretization& time_discretization,
                        Direction& d);

  ///
  /// @brief Gets the maximum primal step size of the fraction-to-boundary-rule.
  /// @return The primal step size of the fraction-to-boundary-rule.
  ///
  double maxPrimalStepSize() const;

  ///
  /// @brief Gets the maximum dual step size of the fraction-to-boundary-rule.
  /// @return The dual step size of the fraction-to-boundary-rule.
  ///
  double maxDualStepSize() const;

  ///
  /// @brief Integrates the solution. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] dual_step_size Dual step size.
  /// @param[in, out] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integrateSolution(const aligned_vector<Robot>& robots,
                         const TimeDiscretization& time_discretization, 
                         const double primal_step_size,
                         const double dual_step_size,
                         Direction& d, Solution& s);

  ///
  /// @brief Integrates the primal solution. 
  /// @param[in, out] robots aligned_vector of Robot for paralle computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] primal_step_size Primal step size.
  /// @param[in] d Direction. 
  /// @param[in, out] s Solution. 
  ///
  void integratePrimalSolution(const aligned_vector<Robot>& robots,
                               const TimeDiscretization& time_discretization, 
                               const double primal_step_size,
                               const Direction& d, Solution& s);

  ///
  /// @brief Resizes the internal data. 
  /// @param[in] time_discretization Time discretization. 
  ///
  void resizeData(const TimeDiscretization& time_discretization);

private:
  int nthreads_;
  aligned_vector<OCPData> ocp_data_;
  IntermediateStage intermediate_stage_;
  ImpactStage impact_stage_;
  TerminalStage terminal_stage_;
  PerformanceIndex performance_index_; 
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;
};

} // namespace robotoc 

#endif // ROBOTOC_DIRECT_MULTIPLE_SHOOTING_HPP_ 