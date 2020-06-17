#ifndef IDOCP_LINE_SEARCH_FILTER_HPP_
#define IDOCP_LINE_SEARCH_FILTER_HPP_

#include <vector>
#include <utility>


namespace idocp {

class LineSearchFilter {
public:
  LineSearchFilter();

  bool isAccepted(const double cost, const double constraint_violation);

  void append(const double cost, const double constraint_violation);

  void clear();
 
private:
  std::vector<std::pair<double, double>> filter_;
  double cost_reduction_rate_, constraints_reduction_rate_;

};

} // namespace idocp


#endif // IDOCP_LINE_SEARCH_FILTER_HPP_