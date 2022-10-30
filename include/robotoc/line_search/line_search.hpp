#ifndef ROBOTOC_LINE_SEARCH_HPP_
#define ROBOTOC_LINE_SEARCH_HPP_

#include <memory>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/core/solution.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/ocp/direct_multiple_shooting.hpp"
#include "robotoc/line_search/line_search_filter.hpp"
#include "robotoc/line_search/line_search_settings.hpp"


namespace robotoc {

///
/// @class LineSearch 
/// @brief Line search for optimal control problems.
///
class LineSearch {
public:
  ///
  /// @brief Construct a line search.
  /// @param[in] ocp Optimal control problem. 
  /// @param[in] settings Line search settings.
  ///
  LineSearch(const OCP& ocp, 
             const LineSearchSettings& settings=LineSearchSettings());

  ///
  /// @brief Default constructor. 
  ///
  LineSearch();

  ///
  /// @brief Default destructor. 
  ///
  ~LineSearch() = default;

  ///
  /// @brief Default copy constructor. 
  ///
  LineSearch(const LineSearch&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  LineSearch& operator=(const LineSearch&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  LineSearch(LineSearch&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  LineSearch& operator=(LineSearch&&) noexcept = default;

  ///
  /// @brief Compute primal step size by fliter line search method. 
  /// @param[in] dms Direct multiple shooting structure.
  /// @param[in, out] robots aligned_vector of Robot for parallel computing.
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] q Initial configuration.
  /// @param[in] v Initial generalized velocity.
  /// @param[in] s Solution. 
  /// @param[in] d Direction. 
  /// @param[in] max_primal_step_size Maximum primal step size. 
  ///
  double computeStepSize(
      const DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
      const TimeDiscretization& time_discretization,
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
      const Direction& d, const double max_primal_step_size);

  ///
  /// @brief Clear the line search filter. 
  ///
  void clearHistory();

  ///
  /// @brief Set line search settings.
  /// @param[in] settings Line search settings. 
  ///
  void set(const LineSearchSettings& settings);

  ///
  /// @brief Resizes the internal data. 
  /// @param[in] time_discretization Time discretization. 
  ///
  void resizeData(const TimeDiscretization& time_discretization);

private:
  LineSearchFilter filter_;
  LineSearchSettings settings_;
  DirectMultipleShooting dms_trial_;
  Solution s_trial_;
  KKTResidual kkt_residual_;

  double lineSearchFilterMethod(
      const DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
      const TimeDiscretization& time_discretization,
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
      const Direction& d, const double max_primal_step_size);

  double meritBacktrackingLineSearch(
      const DirectMultipleShooting& dms, aligned_vector<Robot>& robots,
      const TimeDiscretization& time_discretization,
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
      const Direction& d, const double max_primal_step_size);

  bool armijoCondition(const double merit, const double merit_trial, 
                       const double merit_directional_derivative, 
                       const double step_size) const;

  double penaltyParam(const TimeDiscretization time_discretization, 
                      const Solution& s) const;

};

} // namespace robotoc 

#endif // ROBOTOC_LINE_SEARCH_HPP_ 