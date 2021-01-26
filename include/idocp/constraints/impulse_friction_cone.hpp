#ifndef IDOCP_IMPULSE_FRICTION_CONE_HPP_
#define IDOCP_IMPULSE_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/constraints/impulse_constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"


namespace idocp {

class ImpulseFrictionCone final : public ImpulseConstraintComponentBase {
public:
  ImpulseFrictionCone(const Robot& robot, const double barrier=1.0e-04,
                      const double fraction_to_boundary_rate=0.995);

  ImpulseFrictionCone();

  ~ImpulseFrictionCone();

  // Use default copy constructor.
  ImpulseFrictionCone(const ImpulseFrictionCone&) = default;

  // Use default copy coperator.
  ImpulseFrictionCone& operator=(const ImpulseFrictionCone&) = default;

  // Use default move constructor.
  ImpulseFrictionCone(ImpulseFrictionCone&&) noexcept = default;

  // Use default move assign coperator.
  ImpulseFrictionCone& operator=(ImpulseFrictionCone&&) noexcept = default;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const;

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

  static double frictionConeResidual(const double mu, 
                                     const Eigen::Vector3d& f) {
    assert(mu > 0);
    return (f.coeff(0)*f.coeff(0) + f.coeff(1)*f.coeff(1)
            - mu*mu*f.coeff(2)*f.coeff(2));
  }

private:
  int dimc_;
  double fraction_to_boundary_rate_;

};

} // namespace idocp

#endif // IDOCP_IMPULSE_FRICTION_CONE_HPP_ 