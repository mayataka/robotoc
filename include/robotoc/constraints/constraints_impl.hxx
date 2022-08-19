#ifndef ROBOTOC_CONSTRAINTS_IMPL_HXX_
#define ROBOTOC_CONSTRAINTS_IMPL_HXX_

#include "robotoc/constraints/constraints_impl.hpp"

#include <cassert>
#include <cmath>

namespace robotoc {
namespace constraintsimpl {

template <typename ConstraintComponentBaseTypePtr>
inline void clear(std::vector<ConstraintComponentBaseTypePtr>& constraints) {
  constraints.clear();
}


template <typename ConstraintComponentBaseTypePtr>
inline void createConstraintsData(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    std::vector<ConstraintComponentData>& data) {
  data.clear();
  for (const auto& constraint : constraints) {
    auto component_data = ConstraintComponentData(constraint->dimc(), 
                                                  constraint->barrier());
    constraint->allocateExtraData(component_data);
    data.push_back(component_data);
  }
}


template <typename ConstraintComponentBaseTypePtr, typename ContactStatusType, 
          typename SplitSolutionType>
inline bool isFeasible(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    Robot& robot, const ContactStatusType& contact_status, 
    std::vector<ConstraintComponentData>& data, const SplitSolutionType& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    bool feasible = constraints[i]->isFeasible(robot, contact_status, data[i], s);
    if (!feasible) {
      return false;
    }
  }
  return true;
}


template <typename ConstraintComponentBaseTypePtr, typename ContactStatusType, 
          typename SplitSolutionType>
inline void setSlackAndDual(
   const std::vector<ConstraintComponentBaseTypePtr>& constraints,
   Robot& robot, const ContactStatusType& contact_status, 
   std::vector<ConstraintComponentData>& data, const SplitSolutionType& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->setSlack(robot, contact_status, data[i], s);
    constraints[i]->setSlackAndDualPositive(data[i]);
  }
}


template <typename ConstraintComponentBaseTypePtr, typename ContactStatusType, 
          typename SplitSolutionType>
inline void evalConstraint(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, const ContactStatusType& contact_status, 
    std::vector<ConstraintComponentData>& data, const SplitSolutionType& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->evalConstraint(robot, contact_status, data[i], s);
  }
}


template <typename ConstraintComponentBaseTypePtr, typename ContactStatusType, 
          typename SplitSolutionType, typename SplitKKTResidualType>
inline void linearizeConstraints(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, const ContactStatusType& contact_status, 
    std::vector<ConstraintComponentData>& data, const SplitSolutionType& s, 
    SplitKKTResidualType& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->evalConstraint(robot, contact_status, data[i], s);
    constraints[i]->evalDerivatives(robot, contact_status, data[i], s, 
                                    kkt_residual);
  }
}


template <typename ConstraintComponentBaseTypePtr, typename ContactStatusType,
          typename SplitKKTMatrixType, typename SplitKKTResidualType>
inline void condenseSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    const ContactStatusType& contact_status,
    std::vector<ConstraintComponentData>& data, SplitKKTMatrixType& kkt_matrix, 
    SplitKKTResidualType& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->condenseSlackAndDual(contact_status, data[i], kkt_matrix, 
                                         kkt_residual);
  }
}


template <typename ConstraintComponentBaseTypePtr, typename ContactStatusType,
          typename SplitDirectionType>
inline void expandSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const ContactStatusType& contact_status, 
    std::vector<ConstraintComponentData>& data, const SplitDirectionType& d) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->expandSlackAndDual(contact_status, data[i], d);
  }
}


template <typename ConstraintComponentBaseTypePtr>
inline double maxSlackStepSize(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double min_step_size = 1;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    const double step_size = constraints[i]->maxSlackStepSize(data[i]);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


template <typename ConstraintComponentBaseTypePtr>
inline double maxDualStepSize(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double min_step_size = 1;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    const double step_size = constraints[i]->maxDualStepSize(data[i]);
    if (step_size < min_step_size) {
      min_step_size = step_size;
    }
  }
  return min_step_size;
}


inline void updateSlack(std::vector<ConstraintComponentData>& data, 
                        const double step_size) {
  assert(step_size > 0);
  for (auto& e : data) {
    e.slack.noalias() += step_size * e.dslack;
  }
}


inline void updateDual(std::vector<ConstraintComponentData>& data, 
                       const double step_size) {
  assert(step_size > 0);
  for (auto& e : data) {
    e.dual.noalias() += step_size * e.ddual;
  }
}


template <typename ConstraintComponentBaseTypePtr>
inline void setBarrier(std::vector<ConstraintComponentBaseTypePtr>& constraints, 
                       const double barrier) {
  for (auto& constraint : constraints) {
    constraint->setBarrier(barrier);
  }
}


template <typename ConstraintComponentBaseTypePtr>
inline void setFractionToBoundaryRule(
    std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const double fraction_to_boundary_rule) {
  for (auto& constraint : constraints) {
    constraint->setFractionToBoundaryRule(fraction_to_boundary_rule);
  }
}

} // namespace constraintsimpl
} // namespace robotoc

#endif // ROBOTOC_CONSTRAINTS_IMPL_HXX_ 