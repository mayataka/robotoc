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

using ConstraintComponentBasePtr = std::shared_ptr<ConstraintComponentBase>;
using ImpulseConstraintComponentBasePtr 
  = std::shared_ptr<ImpulseConstraintComponentBase>;

template <typename ConstraintComponentBaseTypePtr>
void clear(std::vector<ConstraintComponentBaseTypePtr>& constraints);

template <typename ConstraintComponentBaseTypePtr>
bool useKinematics(
   const std::vector<ConstraintComponentBaseTypePtr>& constraints);

template <typename ConstraintComponentBaseTypePtr>
void createConstraintsData(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints, 
    std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
bool isFeasible(const std::vector<ConstraintComponentBaseTypePtr>& constraints,
                Robot& robot, std::vector<ConstraintComponentData>& data, 
                const SplitSolutionType& s);

template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
void setSlackAndDual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s);

void augmentDualResidual(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const double dtau, 
    const SplitSolution& s, SplitKKTResidual& kkt_residual);

void augmentDualResidual(
    const std::vector<ImpulseConstraintComponentBasePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTResidual& kkt_residual);

void condenseSlackAndDual(
    const std::vector<ConstraintComponentBasePtr>& constraints, Robot& robot, 
    std::vector<ConstraintComponentData>& data, const double dtau, 
    const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
    SplitKKTResidual& kkt_residual);

void condenseSlackAndDual(
    const std::vector<ImpulseConstraintComponentBasePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
    ImpulseSplitKKTResidual& kkt_residual);

template <typename ConstraintComponentBaseTypePtr, 
          typename SplitSolutionType, typename SplitDirectionType>
void computeSlackAndDualDirection(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s, const SplitDirectionType& d);

template <typename ConstraintComponentBaseTypePtr>
double maxSlackStepSize(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseTypePtr>
double maxDualStepSize(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data);

void updateSlack(std::vector<ConstraintComponentData>& data, 
                 const double step_size);

void updateDual(std::vector<ConstraintComponentData>& data, 
                const double step_size);

template <typename ConstraintComponentBaseTypePtr>
double costSlackBarrier(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseTypePtr>
double costSlackBarrier(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data, 
    const double step_size);

template <typename ConstraintComponentBaseTypePtr, typename SplitSolutionType>
void computePrimalAndDualResidual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    Robot& robot, std::vector<ConstraintComponentData>& data, 
    const SplitSolutionType& s);

template <typename ConstraintComponentBaseTypePtr>
double l1NormPrimalResidual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data);

template <typename ConstraintComponentBaseTypePtr>
double squaredNormPrimalAndDualResidual(
    const std::vector<ConstraintComponentBaseTypePtr>& constraints,
    const std::vector<ConstraintComponentData>& data);


double l1NormPrimalResidual(const std::vector<ConstraintComponentData>& data);

double squaredNormPrimalAndDualResidual(
    const std::vector<ConstraintComponentData>& data);

} // namespace constraintsimpl
} // namespace idocp

#include "idocp/constraints/constraints_impl.hxx"

#endif // IDOCP_CONSTRAINTS_IMPL_HPP_