#ifndef IDOCP_LINE_SEARCH_FILTER_HPP_
#define IDOCP_LINE_SEARCH_FILTER_HPP_

#include <vector>
#include <utility>


namespace idocp {

///
/// @class LineSearchFilter
/// @brief Filter of the line search method to reduce the constraint violation. 
///
class LineSearchFilter {
public:
  LineSearchFilter(const double cost_reduction_rate=0.005, 
                   const double constraints_reduction_rate=0.005);

  ~LineSearchFilter();

  ///
  /// @brief Default copy constructor. 
  ///
  LineSearchFilter(const LineSearchFilter&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  LineSearchFilter& operator=(const LineSearchFilter&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  LineSearchFilter(LineSearchFilter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  LineSearchFilter& operator=(LineSearchFilter&&) noexcept = default;

  bool isAccepted(const double cost, const double constraint_violation);

  void augment(const double cost, const double constraint_violation);

  void clear();

  bool isEmpty() const;
 
private:
  std::vector<std::pair<double, double>> filter_;
  double cost_reduction_rate_, constraints_reduction_rate_;

};

} // namespace idocp


#endif // IDOCP_LINE_SEARCH_FILTER_HPP_