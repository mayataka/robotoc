#ifndef IDOCP_TEST_HELPER_HPP_
#define IDOCP_TEST_HELPER_HPP_

#include <vector>
#include <memory>
#include <cassert>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/discrete_event.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"
#include "idocp/hybrid/parnmpc_discretizer.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {
namespace testhelper {

ContactSequence CreateContactSequence(const Robot& robot, const int N, 
                                      const int max_num_impulse,
                                      const double t0,
                                      const double event_period);

std::shared_ptr<CostFunction> CreateCost(const Robot& robot);

std::shared_ptr<Constraints> CreateConstraints(const Robot& robot);

Solution CreateSolution(const Robot& robot, const int N, const int max_num_impulse=0);

Solution CreateSolution(const Robot& robot, const ContactSequence& contact_sequence, 
                        const double T, const int N, const int max_num_impulse, const double t, 
                        const bool is_parnmpc=false);

Direction CreateDirection(const Robot& robot, const int N, const int max_num_impulse=0);

Direction CreateDirection(const Robot& robot, const ContactSequence& contact_sequence, 
                          const double T, const int N, const int max_num_impulse, const double t,
                          const bool is_parnmpc=false);

KKTMatrix CreateKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence, 
                          const int N, const int max_num_impulse, const bool is_parnmpc=false);

KKTResidual CreateKKTResidual(const Robot& robot, const ContactSequence& contact_sequence, 
                              const int N, const int max_num_impulse, const bool is_parnmpc=false);

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