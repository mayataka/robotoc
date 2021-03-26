#ifndef IDOCP_TEST_HELPER_HPP_
#define IDOCP_TEST_HELPER_HPP_

#include "Eigen/Core"

#include "idocp/hybrid/hybrid_container.hpp"


namespace idocp {
namespace testhelper {

template <typename Type, typename ImpulseType>
bool IsApprox(const hybrid_container<Type, ImpulseType>& rhs, 
              const hybrid_container<Type, ImpulseType>& lhs) {
  assert(rhs.data.size() == lhs.data.size());
  assert(rhs.impulse.size() == lhs.impulse.size());
  assert(rhs.aux.size() == lhs.aux.size());
  assert(rhs.lift.size() == lhs.lift.size());
  const int N = rhs.data.size()-1;
  const int max_num_impulse = rhs.impulse.size();
  for (int i=0; i<=N; ++i) {
    if (!rhs[i].isApprox(lhs[i])) {
      return false;
    } 
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (!rhs.impulse[i].isApprox(lhs.impulse[i])) {
      return false;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (!rhs.aux[i].isApprox(lhs.aux[i])) {
      return false;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (!rhs.lift[i].isApprox(lhs.lift[i])) {
      return false;
    }
  }
  return true;
}


template <typename Type, typename ImpulseType>
bool HasNaN(const hybrid_container<Type, ImpulseType>& obj) {
  const int N = obj.data.size()-1;
  const int max_num_impulse = obj.impulse.size();
  for (int i=0; i<=N; ++i) {
    if (obj[i].hasNaN()) {
      return true;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (obj.impulse[i].hasNaN()) {
      return true;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (obj.aux[i].hasNaN()) {
      return true;
    }
  }
  for (int i=0; i<max_num_impulse; ++i) {
    if (obj.lift[i].hasNaN()) {
      return true;
    }
  }
  return false;
}

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_HPP_