#ifndef ROBOTOC_TEST_HELPER_HPP_
#define ROBOTOC_TEST_HELPER_HPP_

#include "Eigen/Core"


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


template <typename Type, typename ImpulseType, typename SwitchingType>
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