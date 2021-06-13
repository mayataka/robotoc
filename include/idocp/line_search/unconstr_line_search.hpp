#ifndef IDOCP_UNCONSTR_LINE_SEARCH_HPP_
#define IDOCP_UNCONSTR_LINE_SEARCH_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/utils/aligned_vector.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/unconstr/unconstr_ocp.hpp"
#include "idocp/unconstr/unconstr_parnmpc.hpp"
#include "idocp/line_search/line_search_filter.hpp"


namespace idocp {

///
/// @class UnconstrLineSearch 
/// @brief Line search for optimal control problems for unconstrained 
/// rigid-body systems.
///
class UnconstrLineSearch {
public:
  ///
  /// @brief Construct a line search.
  /// @param[in] robot Robot model. 
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  /// @param[in] step_size_reduction_rate Reduction rate of the step size. 
  /// Defalt is 0.75.
  /// @param[in] min_step_size Minimum step size. Default is 0.05.
  ///
  UnconstrLineSearch(const Robot& robot, const double T, const int N, 
                     const int nthreads=1, 
                     const double step_size_reduction_rate=0.75, 
                     const double min_step_size=0.05);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrLineSearch();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrLineSearch();

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrLineSearch(const UnconstrLineSearch&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnconstrLineSearch& operator=(const UnconstrLineSearch&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrLineSearch(UnconstrLineSearch&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrLineSearch& operator=(UnconstrLineSearch&&) noexcept = default;

  ///
  /// @brief Compute primal step size by fliter line search method. 
  /// @param[in, out] ocp optimal control problem.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] t Initial time of the horizon. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in] d Direction. 
  /// @param[in] max_primal_step_size Maximum primal step size. 
  ///
  template <typename UnconstrOCPType>
  double computeStepSize(UnconstrOCPType& ocp, aligned_vector<Robot>& robots, 
                         const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v, const Solution& s, 
                         const Direction& d, 
                         const double max_primal_step_size) {
    assert(max_primal_step_size > 0);
    assert(max_primal_step_size <= 1);
    // If filter is empty, augment the current solution to the filter.
    if (filter_.isEmpty()) {
      computeCostAndViolation(ocp, robots, t, q, v, s);
      filter_.augment(totalCosts(), totalViolations());
    }
    double primal_step_size = max_primal_step_size;
    while (primal_step_size > min_step_size_) {
      computeSolutionTrial(s, d, primal_step_size);
      computeCostAndViolation(ocp, robots, t, q, v, s_trial_, primal_step_size);
      const double total_costs = totalCosts();
      const double total_violations = totalViolations();
      if (filter_.isAccepted(total_costs, total_violations)) {
        filter_.augment(total_costs, total_violations);
        break;
      }
      primal_step_size *= step_size_reduction_rate_;
    }
    if (primal_step_size > min_step_size_) {
      return primal_step_size;
    }
    else {
      return min_step_size_;
    }
  }

  ///
  /// @brief Clear the line search filter. 
  ///
  void clearFilter();

  ///
  /// @brief Clear the line search filter. 
  ///
  bool isFilterEmpty() const;

private:
  LineSearchFilter filter_;
  int N_, nthreads_;
  double T_, dt_, step_size_reduction_rate_, min_step_size_;
  Eigen::VectorXd costs_, violations_;
  Solution s_trial_;
  KKTResidual kkt_residual_;

  void computeCostAndViolation(UnconstrOCP& ocp, aligned_vector<Robot>& robots,
                               const double t, const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, const Solution& s,
                               const double primal_step_size_for_barrier=0);

  void computeCostAndViolation(UnconstrParNMPC& parnmpc, aligned_vector<Robot>& robots,
                               const double t, const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, const Solution& s,
                               const double primal_step_size_for_barrier=0);

  void computeSolutionTrial(const Solution& s, const Direction& d, 
                            const double step_size);

  static void computeSolutionTrial(const SplitSolution& s, 
                                  const SplitDirection& d, 
                                  const double step_size, 
                                  SplitSolution& s_trial) {
    s_trial.q = s.q + step_size * d.dq();
    s_trial.v = s.v + step_size * d.dv();
    s_trial.a = s.a + step_size * d.da();
    s_trial.u = s.u + step_size * d.du;
  }

  void clearCosts() {
    costs_.setZero();
  }

  void clearViolations() {
    violations_.setZero();
  }

  double totalCosts() const {
    return costs_.sum();
  }

  double totalViolations() const {
    return violations_.sum();
  }

};

} // namespace idocp 

#endif // IDOCP_UNCONSTR_LINE_SEARCH_HPP_ 