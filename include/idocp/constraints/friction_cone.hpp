#ifndef IDOCP_FRICTION_CONE_HPP_
#define IDOCP_FRICTION_CONE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


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

  bool useKinematics() const override;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const double dtau, const SplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const double dtau, const SplitSolution& s,
                           KKTResidual& kkt_residual) const override;

  void augmentDualResidual(const Robot& robot, ConstraintComponentData& data, 
                           const double dtau, const Eigen::VectorXd& u,
                           Eigen::VectorXd& lu) const override {}

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const SplitSolution& s,
                            KKTMatrix& kkt_matrix,
                            KKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const Eigen::VectorXd& u,
                            Eigen::MatrixXd& Quu, 
                            Eigen::VectorXd& lu) const override {}

  void computeSlackAndDualDirection(Robot& robot, ConstraintComponentData& data, 
                                    const double dtau, const SplitSolution& s,
                                    const SplitDirection& d) const override; 

  double residualL1Nrom(Robot& robot, ConstraintComponentData& data, 
                        const double dtau, 
                        const SplitSolution& s) const override;

  double squaredKKTErrorNorm(Robot& robot, ConstraintComponentData& data, 
                             const double dtau, 
                             const SplitSolution& s) const override;
  
  int dimc() const override;

private:
  int dimc_;
  double mu_;

};

} // namespace idocp

#endif // IDOCP_FRICTION_CONE_HPP_ 