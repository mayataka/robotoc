#ifndef ROBOTOC_TEST_HELPER_HPP_
#define ROBOTOC_TEST_HELPER_HPP_

#include "Eigen/Core"

#include "robotoc/utils/aligned_vector.hpp"


namespace robotoc {
namespace testhelper {

template <typename Type>
bool IsApprox(const aligned_vector<Type>& rhs, 
              const aligned_vector<Type>& lhs) {
  assert(rhs.size() == lhs.size());
  for (int i=0; i<rhs.size(); ++i) {
    std::cout << "i = " << i << std::endl;
    if (!rhs[i].isApprox(lhs[i])) {
      return false;
    } 
  }
  return true;
}


template <typename Type>
bool HasNaN(const aligned_vector<Type>& obj) {
  for (int i=0; i<obj.size(); ++i) {
    std::cout << "i = " << i << std::endl;
    if (obj[i].hasNaN()) {
      return true;
    }
  }
  return false;
}

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_HPP_