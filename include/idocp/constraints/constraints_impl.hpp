#ifndef IDOCP_CONSTRAINTS_IMPL_HPP_
#define IDOCP_CONSTRAINTS_IMPL_HPP_

#include <vector>
#include <memory>

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/impulse_constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/constraints/constraints_data.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {
namespace constraintsimpl {

template <typename ConstraintComponentBaseType>
void clear(std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints);

template <typename ConstraintComponentBaseType>
bool useKinematics(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints);

template <typename ConstraintComponentBaseType, typename SplitSolutionType>
bool isFeasible(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const SplitSolutionType& s);

void setSlackAndDual(
   const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const double dtau, const SplitSolution& s);

void setSlackAndDual(
   const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const ImpulseSplitSolution& s);

void augmentDualResidual(
   const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const double dtau, const SplitSolution& s, SplitKKTResidual& kkt_residual);

void augmentDualResidual(
   const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual);

void condenseSlackAndDual(
   const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const double dtau, const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
   SplitKKTResidual& kkt_residual);

void condenseSlackAndDual(
   const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
   ImpulseSplitKKTResidual& kkt_residual);

void computeSlackAndDualDirection(
   const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const double dtau, const SplitSolution& s, const SplitDirection& d);

void computeSlackAndDualDirection(
   const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const ImpulseSplitSolution& s, const ImpulseSplitDirection& d);

template <typename ConstraintComponentBaseType>
double maxSlackStepSize(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   const std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseType>
double maxDualStepSize(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   const std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseType>
void updateSlack(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   std::vector<ConstraintComponentData>& data, const double step_size);

template <typename ConstraintComponentBaseType>
void updateDual(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   std::vector<ConstraintComponentData>& data, const double step_size);

template <typename ConstraintComponentBaseType>
double costSlackBarrier(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   const std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseType>
double costSlackBarrier(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   const std::vector<ConstraintComponentData>& data, 
   const double step_size);

void computePrimalAndDualResidual(
   const std::vector<std::shared_ptr<ConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const double dtau, const SplitSolution& s);

void computePrimalAndDualResidual(
   const std::vector<std::shared_ptr<ImpulseConstraintComponentBase>>& constraints,
   Robot& robot, std::vector<ConstraintComponentData>& data, 
   const ImpulseSplitSolution& s);

template <typename ConstraintComponentBaseType>
double l1NormPrimalResidual(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   const std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseType>
double squaredNormPrimalAndDualResidual(
   const std::vector<std::shared_ptr<ConstraintComponentBaseType>>& constraints,
   const std::vector<ConstraintComponentData>& data);

} // namespace constraintsimpl
} // namespace idocp

#include "idocp/constraints/constraints_impl.hxx"

#endif // IDOCP_CONSTRAINTS_IMPL_HPP_