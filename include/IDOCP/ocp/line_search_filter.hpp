#ifndef IDOCP_LINE_SEARCH_FILTER_HPP_
#define IDOCP_LINE_SEARCH_FILTER_HPP_

#include <vector>
#include <utility>


namespace idocp {

class LineSearchFilter {
public:
  LineSearchFilter();

  bool isAccepted(const double cost, const double constraint_residual);

  void append(const double cost, const double constraint_residual);

  void clear();
 
private:
  std::vector<std::pair<double, double>> filter_;

};

} // namespace idocp


#endif // IDOCP_LINE_SEARCH_FILTER_HPP_