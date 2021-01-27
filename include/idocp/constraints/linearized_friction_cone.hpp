#ifndef IDOCP_LINEARIZED_FRICTION_CONE_HPP_
#define IDOCP_LINEARIZED_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class LinearizedFrictionCone final : public ConstraintComponentBase {
public:
  LinearizedFrictionCone(const Robot& robot, const double mu, 
                         const double barrier=1.0e-04,
                         const double fraction_to_boundary_rate=0.995);

  LinearizedFrictionCone();

  ~LinearizedFrictionCone();

  // Use default copy constructor.
  LinearizedFrictionCone(const LinearizedFrictionCone&) = default;

  // Use default copy coperator.
  LinearizedFrictionCone& operator=(const LinearizedFrictionCone&) = default;

  // Use default move constructor.
  LinearizedFrictionCone(LinearizedFrictionCone&&) noexcept = default;

  // Use default move assign coperator.
  LinearizedFrictionCone& operator=(LinearizedFrictionCone&&) noexcept = default;

  void setFrictionCoefficient(const double mu);

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const SplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const double dtau, const SplitSolution& s,
                           SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const SplitSolution& s,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void computeSlackAndDualDirection(Robot& robot, ConstraintComponentData& data, 
                                    const SplitSolution& s,
                                    const SplitDirection& d) const override; 

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const SplitSolution& s) const override;

  int dimc() const override;

  static Eigen::Vector2d frictionConeResidual(const double mu, 
                                              const Eigen::Vector3d& f) {
    assert(mu > 0);
    Eigen::Vector2d res;
    res.coeffRef(0) = f.coeff(0) - mu * f.coeff(2) / std::sqrt(2);
    res.coeffRef(1) = f.coeff(1) - mu * f.coeff(2) / std::sqrt(2);
    return res;
  }

  static double normalForceResidual(const Eigen::Vector3d& f) {
    return (- f.coeff(2));
  }

private:
  int dimc_;
  double mu_;
  Eigen::Matrix3d Jac;

};

} // namespace idocp

#endif // IDOCP_LINEARIZED_FRICTION_CONE_HPP_ 