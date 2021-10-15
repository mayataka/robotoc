#ifndef IDOCP_LINE_SEARCH_FILTER_HPP_
#define IDOCP_LINE_SEARCH_FILTER_HPP_

#include <set>
#include <utility>


namespace idocp {

///
/// @class LineSearchFilter
/// @brief Filter of the line search. 
///
class LineSearchFilter {
public:
  ///
  /// @brief Construct a line search filter.
  /// @param[in] cost_reduction_rate Reduction rate of the cost. 
  /// @param[in] constraints_reduction_rate Reduction rate of the constraints. 
  ///
  LineSearchFilter(const double cost_reduction_rate=0.005, 
                   const double constraints_reduction_rate=0.005);

  ///
  /// @brief Destructor. 
  ///
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

  ///
  /// @brief Check wheather the pair of cost and constraint violation is 
  /// accepted by the filter or not. 
  /// @param[in] cost Value of the cost.
  /// @param[in] constraint_violation Value of the constraint violation. 
  /// @return true if the pair of cost and constraint violation is accepted.
  /// false if not. 
  ///
  bool isAccepted(const double cost, const double constraint_violation);

  ///
  /// @brief Augment the pair of cost and constraint violation to the filter. 
  /// @param[in] cost Value of the cost.
  /// @param[in] constraint_violation Value of the constraint violation. 
  ///
  void augment(const double cost, const double constraint_violation);

  ///
  /// @brief Clears the filter. 
  ///
  void clear();

  ///
  /// @brief Checks wheather the filter is empty or not. 
  /// @return true if the filter is empty. false if not.
  ///
  bool isEmpty() const;

private:
  std::set<std::pair<double, double>> filter_;
  double cost_reduction_rate_, constraints_reduction_rate_;

};

} // namespace idocp


#endif // IDOCP_LINE_SEARCH_FILTER_HPP_