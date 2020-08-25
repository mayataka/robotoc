#ifndef IDOCP_JOINT_TORQUES_UPPER_LIMIT_HPP_
#define IDOCP_JOINT_TORQUES_UPPER_LIMIT_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/kkt_matrix.hpp"


namespace idocp {

class JointTorquesUpperLimit final : public ConstraintComponentBase {
public:
  JointTorquesUpperLimit(const Robot& robot, const double barrier=1.0e-04,
                          const double fraction_to_boundary_rate=0.995);

  JointTorquesUpperLimit();

  ~JointTorquesUpperLimit();

  // Use default copy constructor.
  JointTorquesUpperLimit(const JointTorquesUpperLimit&) = default;

  // Use default copy coperator.
  JointTorquesUpperLimit& operator=(const JointTorquesUpperLimit&) = default;

  // Use default move constructor.
  JointTorquesUpperLimit(JointTorquesUpperLimit&&) noexcept = default;

  // Use default move assign coperator.
  JointTorquesUpperLimit& operator=(JointTorquesUpperLimit&&) noexcept 
      = default;

  bool isFeasible(const Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(const Robot& robot, ConstraintComponentData& data, 
                       const double dtau, const SplitSolution& s) const override;

  void augmentDualResidual(const Robot& robot, ConstraintComponentData& data, 
                           const double dtau, 
                           KKTResidual& kkt_residual) const override {}

  void augmentDualResidual(const Robot& robot, ConstraintComponentData& data, 
                           const double dtau, 
                           Eigen::VectorXd& lu) const override;

  void condenseSlackAndDual(const Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const SplitSolution& s,
                            KKTMatrix& kkt_matrix,
                            KKTResidual& kkt_residual) const override {}

  void condenseSlackAndDual(const Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const Eigen::VectorXd& u,
                            Eigen::MatrixXd& Quu, 
                            Eigen::VectorXd& lu) const override;

  void computeSlackAndDualDirection(const Robot& robot, 
                                    ConstraintComponentData& data, 
                                    const double dtau, 
                                    const SplitDirection& d) const override; 

  double residualL1Nrom(const Robot& robot, ConstraintComponentData& data, 
                        const double dtau, 
                        const SplitSolution& s) const override;

  double squaredKKTErrorNorm(const Robot& robot, ConstraintComponentData& data, 
                             const double dtau, 
                             const SplitSolution& s) const override;

  int dimc() const override;

private:
  int dimc_, dim_passive_;
  Eigen::VectorXd umax_;

};

} // namespace idocp

#endif // IDOCP_JOINT_TORQUES_UPPER_LIMIT_HPP_