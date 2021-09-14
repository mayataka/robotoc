#ifndef IDOCP_CONSTRAINTS_IMPL_HXX_
#define IDOCP_CONSTRAINTS_IMPL_HXX_

#include "idocp/constraints/constraints_impl.hpp"

#include <cassert>
#include <cmath>

namespace idocp {
namespace constraintsimpl {

template <typename ConstraintComponentBaseTypePtr>
inline void clear(std::vector<ConstraintComponentBaseTypePtr>& constraints) {
  constraints.clear();
}


template <typename ConstraintComponentBaseTypePtr>
inline bool useKinematics(
   const std::vector<ConstraintComponentBaseTypePtr>& constraints) {
  for (const auto& constraint : constraints) {
    if (constraint->useKinematics()) {
      return true;
    }
  }
  return false;
}


template <typename ConstraintComponentBaseTypePtr>
inline void createConstraintsData(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    std::vector<ConstraintComponentData>& data) {
  data.clear();
  for (const auto& constraint : constraints) {
    auto component_data = ConstraintComponentData(
        constraint->dimc(), constraint->barrierParameter());
    constraint->allocateExtraData(component_data);
    data.push_back(component_data);
  }
}


template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
inline bool isFeasible(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    bool feasible = constraints[i]->isFeasible(robot, data[i], s);
    if (!feasible) {
      return false;
    }
  }
  return true;
}


template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
inline void setSlackAndDual(
   const std::vector<ConstraintComponentBaseTypePtr>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const SplitSolutionType& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->setSlack(robot, data[i], s);
    constraints[i]->setSlackAndDualPositive(data[i]);
  }
}


template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
inline void evalConstraint(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->evalConstraint(robot, data[i], s);
  }
}


inline void linearizeConstraints(
   const std::vector<ConstraintComponentBasePtr>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, const double dt,
   const SplitSolution& s, SplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->evalConstraint(robot, data[i], s);
    constraints[i]->evalDerivatives(robot, data[i], dt, s, kkt_residual);
    kkt_residual.h += data[i].dual.dot(data[i].residual);
  }
}


inline void linearizeConstraints(
   const std::vector<ImpulseConstraintComponentBasePtr>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->evalConstraint(robot, data[i], s);
    constraints[i]->evalDerivatives(robot, data[i], s, kkt_residual);
  }
}


inline void condenseSlackAndDual(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const double dt, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->evalConstraint(robot, data[i], s);
    constraints[i]->evalDerivatives(robot, data[i], dt, s, kkt_residual);
    constraints[i]->condenseSlackAndDual(robot, data[i], dt, s, kkt_matrix, 
                                         kkt_residual);
    kkt_residual.h += data[i].dual.dot(data[i].residual);
    const double r1 = (data[i].residual.array()*data[i].cond.array()).sum();
    const double r2 = (data[i].residual.array()*data[i].cmpl.array()/data[i].slack.array()).sum();
    const double rr = r1 - r2; 
    kkt_residual.h += rr;
    kkt_matrix.Qtt += std::abs(rr) / dt;
  }
}


inline void condenseSlackAndDual(
    const std::vector<ImpulseConstraintComponentBasePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->evalConstraint(robot, data[i], s);
    constraints[i]->evalDerivatives(robot, data[i], s, kkt_residual);
    constraints[i]->condenseSlackAndDual(robot, data[i], s, kkt_matrix, 
                                         kkt_residual);
  }
}


template <typename ConstraintComponentBaseTypePtr, 
          typename SplitSolutionType, typename SplitDirectionType>
inline void expandSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    std::vector<ConstraintComponentData>& data, const double dt, 
    const double dts, const SplitSolutionType& s, const SplitDirectionType& d) {
  for (auto& e : data) {
    e.residual.array() *= (1+(dts/dt));
    e.cmpl.array() *= (1+(dts/dt));
  }
  expandSlackAndDual(constraints, data, s, d);
}


template <typename ConstraintComponentBaseTypePtr, 
          typename SplitSolutionType, typename SplitDirectionType>
inline void expandSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s, const SplitDirectionType& d) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->expandSlackAndDual(data[i], s, d);
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
} // namespace idocp

#endif // IDOCP_CONSTRAINTS_IMPL_HXX_ 