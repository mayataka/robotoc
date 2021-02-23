#ifndef IDOCP_FRICTION_CONE_HPP_
#define IDOCP_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class FrictionCone final : public ConstraintComponentBase {
public:
  FrictionCone(const Robot& robot, const double mu, 
               const double barrier=1.0e-04,
               const double fraction_to_boundary_rate=0.995);

  FrictionCone();

  ~FrictionCone();

  // Use default copy constructor.
  FrictionCone(const FrictionCone&) = default;

  // Use default copy coperator.
  FrictionCone& operator=(const FrictionCone&) = default;

  // Use default move constructor.
  FrictionCone(FrictionCone&&) noexcept = default;

  // Use default move assign coperator.
  FrictionCone& operator=(FrictionCone&&) noexcept = default;

  void setFrictionCoefficient(const double mu);

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const SplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const double dt, const SplitSolution& s,
                           SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const double dt, const SplitSolution& s,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void computeSlackAndDualDirection(Robot& robot, ConstraintComponentData& data, 
                                    const SplitSolution& s,
                                    const SplitDirection& d) const override; 

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const SplitSolution& s) const override;
  
  int dimc() const override;

  static double frictionConeResidual(const double mu, 
                                     const Eigen::Vector3d& f) {
    assert(mu > 0);
    return (f.coeff(0)*f.coeff(0) + f.coeff(1)*f.coeff(1)
            - mu*mu*f.coeff(2)*f.coeff(2));
  }

  static double normalForceResidual(const Eigen::Vector3d& f) {
    return (- f.coeff(2));
  }

private:
  int dimc_;
  double mu_;

};

} // namespace idocp

#endif // IDOCP_FRICTION_CONE_HPP_ 