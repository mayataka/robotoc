#ifndef ROBOTOC_LINE_SEARCH_FILTER_HPP_
#define ROBOTOC_LINE_SEARCH_FILTER_HPP_

#include <vector>
#include <utility>


namespace robotoc {

///
/// @class LineSearchFilter
/// @brief Filter of the line search. 
///
class LineSearchFilter {
public:
  ///
  /// @brief Construct a line search filter.
  /// @param[in] beta A parameter for filter acceptation. 
  ///
  LineSearchFilter(const double beta=1.0e-05);

  ///
  /// @brief Default destructor. 
  ///
  ~LineSearchFilter() = default;

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
  bool isAccepted(const double cost, const double constraint_violation) const;

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
  std::vector<std::pair<double, double>> filter_;
  double beta_;

};

} // namespace robotoc


#endif // ROBOTOC_LINE_SEARCH_FILTER_HPP_