#ifndef IDOCP_IMPULSE_NORMAL_FORCE_HPP_
#define IDOCP_IMPULSE_NORMAL_FORCE_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/constraints/impulse_constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

class ImpulseNormalForce final : public ImpulseConstraintComponentBase {
public:
  ImpulseNormalForce(const Robot& robot, const double barrier=1.0e-04,
                    const double fraction_to_boundary_rate=0.995);

  ImpulseNormalForce();

  ~ImpulseNormalForce();

  // Use default copy constructor.
  ImpulseNormalForce(const ImpulseNormalForce&) = default;

  // Use default copy coperator.
  ImpulseNormalForce& operator=(const ImpulseNormalForce&) = default;

  // Use default move constructor.
  ImpulseNormalForce(ImpulseNormalForce&&) noexcept = default;

  // Use default move assign coperator.
  ImpulseNormalForce& operator=(ImpulseNormalForce&&) noexcept = default;

  KinematicsLevel kinematicsLevel() const override;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const ImpulseSplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const ImpulseSplitSolution& s) const override;

  void augmentDualResidual(
      Robot& robot, ConstraintComponentData& data, 
      const ImpulseSplitSolution& s, 
      ImpulseSplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(
      Robot& robot, ConstraintComponentData& data, 
      const ImpulseSplitSolution& s, ImpulseSplitKKTMatrix& kkt_matrix, 
      ImpulseSplitKKTResidual& kkt_residual) const override;

  void computeSlackAndDualDirection(
      Robot& robot, ConstraintComponentData& data, 
      const ImpulseSplitSolution& s, 
      const ImpulseSplitDirection& d) const override; 

  void computePrimalAndDualResidual(
      Robot& robot, ConstraintComponentData& data, 
      const ImpulseSplitSolution& s) const override;

  int dimc() const override;

private:
  int dimc_;

};

} // namespace idocp

#endif // IDOCP_IMPULSE_NORMAL_FORCE_HPP_ 