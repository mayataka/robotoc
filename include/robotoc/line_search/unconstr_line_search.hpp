#ifndef ROBOTOC_UNCONSTR_LINE_SEARCH_HPP_
#define ROBOTOC_UNCONSTR_LINE_SEARCH_HPP_

#include <vector>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/core/solution.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/unconstr/unconstr_direct_multiple_shooting.hpp"
#include "robotoc/parnmpc/unconstr_backward_correction.hpp"
#include "robotoc/line_search/line_search_filter.hpp"


namespace robotoc {

///
/// @class UnconstrLineSearch 
/// @brief Line search for optimal control problems for unconstrained 
/// rigid-body systems.
///
class UnconstrLineSearch {
public:
  ///
  /// @brief Construct a line search.
  /// @param[in] ocp Optimal control problem. 
  /// @param[in] step_size_reduction_rate Reduction rate of the step size. 
  /// Defalt is 0.75.
  /// @param[in] min_step_size Minimum step size. Default is 0.05.
  ///
  UnconstrLineSearch(const OCP& ocp, const double step_size_reduction_rate=0.75, 
                     const double min_step_size=0.05);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrLineSearch();

  ///
  /// @brief Default destructor. 
  ///
  ~UnconstrLineSearch() = default;

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
  /// @param[in, out] dms Direct multiple shooting method.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in] d Direction. 
  /// @param[in] max_primal_step_size Maximum primal step size. 
  ///
  double computeStepSize(UnconstrDirectMultipleShooting& dms, 
                         aligned_vector<Robot>& robots, 
                         const std::vector<GridInfo>& time_discretization, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         const Solution& s, const Direction& d, 
                         const double max_primal_step_size);

  ///
  /// @brief Compute primal step size by fliter line search method. 
  /// @param[in, out] backward_correction Backward correction method.
  /// @param[in] robots aligned_vector of Robot.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in] d Direction. 
  /// @param[in] max_primal_step_size Maximum primal step size. 
  ///
  double computeStepSize(UnconstrBackwardCorrection& backward_correction, 
                         aligned_vector<Robot>& robots, 
                         const std::vector<GridInfo>& time_discretization, 
                         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                         const Solution& s, const Direction& d, 
                         const double max_primal_step_size);

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
  double step_size_reduction_rate_, min_step_size_;
  Solution s_trial_;
  KKTResidual kkt_residual_;

  void computeSolutionTrial(const Solution& s, const Direction& d, 
                            const double step_size);

  static void computeSolutionTrial(const SplitSolution& s, 
                                  const SplitDirection& d, 
                                  const double step_size, 
                                  SplitSolution& s_trial);

};

} // namespace robotoc 

#endif // ROBOTOC_UNCONSTR_LINE_SEARCH_HPP_ 