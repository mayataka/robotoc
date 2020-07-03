#ifndef IDOCP_LINE_SEARCH_FILTER_HPP_
#define IDOCP_LINE_SEARCH_FILTER_HPP_

#include <vector>
#include <utility>


namespace idocp {

class LineSearchFilter {
public:
  LineSearchFilter();

  ~LineSearchFilter();

  // Use default copy constructor.
  LineSearchFilter(const LineSearchFilter&) = default;

  // Use default copy operator.
  LineSearchFilter& operator=(const LineSearchFilter&) = default;

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