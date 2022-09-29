#ifndef ROBOTOC_UTILS_NUMERICS_HPP_
#define ROBOTOC_UTILS_NUMERICS_HPP_

#include <cmath>
#include <limits>

namespace robotoc {
namespace numerics {

bool isApprox(const double x, const double y, 
              const double eps=std::sqrt(std::numeric_limits<double>::epsilon())) {
  return (std::fabs(x-y) < eps);
}
  
} // namespace numerics
} // namespace robotoc 
;
#endif // ROBOTOC_UTILS_NUMERICS_HPP_