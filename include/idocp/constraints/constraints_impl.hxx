#ifndef IDOCP_CONSTRAINTS_IMPL_HXX_
#define IDOCP_CONSTRAINTS_IMPL_HXX_

#include "idocp/constraints/constraints_impl.hpp"

#include <cassert>

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
    auto component_data = ConstraintComponentData(constraint->dimc());
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
    constraints[i]->setSlackAndDual(robot, data[i], s);
  }
}


inline void augmentDualResidual(
   const std::vector<ConstraintComponentBasePtr>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const double dtau, const SplitSolution& s, SplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->augmentDualResidual(robot, data[i], dtau, s, kkt_residual);
  }
}


inline void augmentDualResidual(
   const std::vector<ImpulseConstraintComponentBasePtr>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->augmentDualResidual(robot, data[i], s, kkt_residual);
  }
}


inline void condenseSlackAndDual(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->condenseSlackAndDual(robot, data[i], dtau, s, kkt_matrix, 
                                         kkt_residual);
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
    constraints[i]->condenseSlackAndDual(robot, data[i], s, kkt_matrix, 
                                         kkt_residual);
  }
}


template <typename ConstraintComponentBaseTypePtr, 
          typename SplitSolutionType, typename SplitDirectionType>
inline void computeSlackAndDualDirection(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s, const SplitDirectionType& d) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computeSlackAndDualDirection(robot, data[i], s, d);
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
inline double costSlackBarrier(
   const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double cost = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    cost += constraints[i]->costSlackBarrier(data[i]);
  }
  return cost;
}


template <typename ConstraintComponentBaseTypePtr>
inline double costSlackBarrier(
   const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data, const double step_size) {
  assert(constraints.size() == data.size());
  double cost = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    cost += constraints[i]->costSlackBarrier(data[i], step_size);
  }
  return cost;
}


template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
inline void computePrimalAndDualResidual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computePrimalAndDualResidual(robot, data[i], s);
  }
}


template <typename ConstraintComponentBaseTypePtr>
inline double l1NormPrimalResidual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double l1_norm = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    l1_norm += constraints[i]->l1NormPrimalResidual(data[i]);
  }
  return l1_norm;
}


template <typename ConstraintComponentBaseTypePtr>
inline double squaredNormPrimalAndDualResidual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double squared_norm = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    squared_norm += data[i].residual.squaredNorm() + data[i].duality.squaredNorm();
  }
  return squared_norm;
}


inline double l1NormPrimalResidual(
    const std::vector<ConstraintComponentData>& data) {
  double l1_norm = 0;
  for (const auto& e : data) {
    l1_norm += e.residual.lpNorm<1>();
  }
  return l1_norm;
}


inline double squaredNormPrimalAndDualResidual(
    const std::vector<ConstraintComponentData>& data) {
  double squared_norm = 0;
  for (const auto& e : data) {
    squared_norm += e.residual.squaredNorm() + e.duality.squaredNorm();
  }
  return squared_norm;
}


template <typename ConstraintComponentBaseTypePtr>
inline void setBarrier(std::vector<ConstraintComponentBaseTypePtr>& constraints, 
                       const double barrier) {
  for (auto& constraint : constraints) {
    constraint->setBarrier(barrier);
  }
}


template <typename ConstraintComponentBaseTypePtr>
inline void setFractionToBoundaryRate(
    std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const double fraction_to_boundary_rate) {
  for (auto& constraint : constraints) {
    constraint->setFractionToBoundaryRate(fraction_to_boundary_rate);
  }
}

} // namespace constraintsimpl
} // namespace idocp

#endif // IDOCP_CONSTRAINTS_IMPL_HXX_ 