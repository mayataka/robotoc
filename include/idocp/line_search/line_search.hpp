#ifndef IDOCP_LINE_SEARCH_HPP_
#define IDOCP_LINE_SEARCH_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/cost/impulse_cost_function.hpp"
#include "idocp/ocp/contact_sequence.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/impulse_split_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/riccati_solution.hpp"
#include "idocp/ocp/line_search_filter.hpp"
#include "idocp/ocp/split_temporary_solution.hpp"
#include "idocp/hybrid/hybrid_container.hpp"


namespace idocp {

///
/// @class LineSearch 
/// @brief Line search for optimal control problems.
///
class LineSearch {
public:
  LineSearch(const int N, const int max_num_impulse=0, const int num_proc=1);

  ///
  /// @brief Default constructor. 
  ///
  LineSearch();

  ///
  /// @brief Destructor. 
  ///
  ~LineSearch();

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

  std::pair<double, double> computeStepSize(
      HybridOCP& split_ocps, std::vector<Robot>& robot,
      const ContactSequence& contact_sequence, const double t, 
      const Eigen::VectorXd& q, const Eigen::VectorXd& v, const Solution& s, 
      const Diretion& d);

  ///
  /// @brief Clear the line search filter. 
  ///
  void clearFilter();

  bool isFilterEmpty() const;

private:
  LineSearchFilter filter_;
  int N_, max_num_impulse_, num_proc_;
  double step_size_reduction_rate_, min_step_size_;
  Solution s_tmp_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_, costs_, violations_;

};

} // namespace idocp 


#endif // IDOCP_LINE_SEARCH_HPP_ 