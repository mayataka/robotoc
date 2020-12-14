#ifndef IDOCP_JOINT_ACCELERATION_UPPER_LIMIT_HPP_
#define IDOCP_JOINT_ACCELERATION_UPPER_LIMIT_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class JointAccelerationUpperLimit final : public ConstraintComponentBase {
public:
  JointAccelerationUpperLimit(const Robot& robot, const Eigen::VectorXd& amax,
                              const double barrier=1.0e-04,
                              const double fraction_to_boundary_rate=0.995);

  JointAccelerationUpperLimit();

  ~JointAccelerationUpperLimit();

  // Use default copy constructor.
  JointAccelerationUpperLimit(const JointAccelerationUpperLimit&) = default;

  // Use default copy coperator.
  JointAccelerationUpperLimit& operator=(const JointAccelerationUpperLimit&) = default;

  // Use default move constructor.
  JointAccelerationUpperLimit(JointAccelerationUpperLimit&&) noexcept = default;

  // Use default move assign coperator.
  JointAccelerationUpperLimit& operator=(JointAccelerationUpperLimit&&) noexcept 
      = default;

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const double dtau, 
                       const SplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const double dtau, const SplitSolution& s, 
                           SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const double dtau, const SplitSolution& s,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void computeSlackAndDualDirection(Robot& robot, ConstraintComponentData& data, 
                                    const double dtau, const SplitSolution& s,
                                    const SplitDirection& d) const override; 

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const double dtau, 
                                    const SplitSolution& s) const override;
  
  int dimc() const override;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

private:
  int dimc_, dim_passive_;
  Eigen::VectorXd amax_;

};

} // namespace idocp

#endif // IDOCP_JOINT_ACCELERATION_UPPER_LIMIT_HPP_