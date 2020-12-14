#ifndef IDOCP_CONSTRAINTS_IMPL_HXX_
#define IDOCP_CONSTRAINTS_IMPL_HXX_

#include "idocp/constraints/constraints_impl.hpp"

#include <cassert>

namespace idocp {
namespace constraintsimpl {

template <typename ConstraintComponentBaseType>
inline void clear(
    std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints) {
  constraints.clear();
}


template <typename ConstraintComponentBaseType>
inline bool useKinematics(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints) {
  for (const auto constraint : constraints) {
    if (constraint->useKinematics()) {
      return true;
    }
  }
  return false;
}


template <typename ConstraintComponentBaseType, typename SplitSolutionType>
inline bool isFeasible(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
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


inline void setSlackAndDual(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->setSlackAndDual(robot, data[i], dtau, s);
  }
}


inline void setSlackAndDual(
    const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const ImpulseSplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->setSlackAndDual(robot, data[i], s);
  }
}


inline void augmentDualResidual(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
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
    const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
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
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
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
    const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
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


template <typename ConstraintComponentBaseType>
inline double maxSlackStepSize(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
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


template <typename ConstraintComponentBaseType>
inline double maxDualStepSize(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
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


template <typename ConstraintComponentBaseType>
inline void updateSlack(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
    std::vector<ConstraintComponentData>& data, const double step_size) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->updateSlack(data[i], step_size);
  }
}


template <typename ConstraintComponentBaseType>
inline void updateDual(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
    std::vector<ConstraintComponentData>& data, const double step_size) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->updateDual(data[i], step_size);
  }
}



template <typename ConstraintComponentBaseType>
inline double costSlackBarrier(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
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


template <typename ConstraintComponentBaseType>
inline double costSlackBarrier(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
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


inline void computePrimalAndDualResidual(
    const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const double dtau, const SplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computePrimalAndDualResidual(robot, data[i], dtau, s);
  }
}


inline void computePrimalAndDualResidual(
    const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const ImpulseSplitSolution& s) {
  assert(constraints.size() == data.size());
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    constraints[i]->computePrimalAndDualResidual(robot, data[i], s);
  }
}


template <typename ConstraintComponentBaseType>
inline double Constraints::l1NormPrimalResidual(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
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


template <typename ConstraintComponentBaseType>
inline double Constraints::squaredNormPrimalAndDualResidual(
    const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
    const std::vector<ConstraintComponentData>& data) {
  assert(constraints.size() == data.size());
  double squared_norm = 0;
  for (int i=0; i<constraints.size(); ++i) {
    assert(data[i].dimc() == constraints[i]->dimc());
    assert(data[i].checkDimensionalConsistency());
    squared_norm += constraints[i]->squaredNormPrimalAndDualResidual(data[i]);
  }
  return squared_norm;
}

} // namespace constraintsimpl
} // namespace idocp

#endif // IDOCP_CONSTRAINTS_IMPL_HXX_ 