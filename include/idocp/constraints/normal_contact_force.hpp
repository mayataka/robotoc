#ifndef IDOCP_NORMAL_CONTACT_FORCE_HPP_
#define IDOCP_NORMAL_CONTACT_FORCE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class NormalContactForce final : public ConstraintComponentBase {
public:
  NormalContactForce(const Robot& robot, const double barrier=1.0e-04, 
                     const double fraction_to_boundary_rate=0.995);

  NormalContactForce();

  ~NormalContactForce();

  // Use default copy constructor.
  NormalContactForce(const NormalContactForce&) = default;

  // Use default copy coperator.
  NormalContactForce& operator=(const NormalContactForce&) = default;

  // Use default move constructor.
  NormalContactForce(NormalContactForce&&) noexcept = default;

  // Use default move assign coperator.
  NormalContactForce& operator=(NormalContactForce&&) noexcept = default;

  bool useKinematics() const override;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const double dtau, const SplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const double dtau, 
                           KKTResidual& kkt_residual) const override;

  void augmentDualResidual(const Robot& robot, ConstraintComponentData& data, 
                           const double dtau, 
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
                                    const double dtau, 
                                    const SplitDirection& d) const override; 

  double residualL1Nrom(Robot& robot, ConstraintComponentData& data, 
                        const double dtau, 
                        const SplitSolution& s) const override;

  double squaredKKTErrorNorm(Robot& robot, ConstraintComponentData& data, 
                             const double dtau, 
                             const SplitSolution& s) const override;
  
  int dimc() const override;

private:
  int contact_index_, dimc_;

};

} // namespace idocp

#endif // IDOCP_NORMAL_CONTACT_FORCE_HPP_