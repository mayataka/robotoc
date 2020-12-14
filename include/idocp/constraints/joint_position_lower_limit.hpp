#ifndef IDOCP_JOINT_POSITION_LOWER_LIMIT_HPP_
#define IDOCP_JOINT_POSITION_LOWER_LIMIT_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

class JointPositionLowerLimit final : public ConstraintComponentBase {
public:
  JointPositionLowerLimit(const Robot& robot, const double barrier=1.0e-04,
                          const double fraction_to_boundary_rate=0.995);

  JointPositionLowerLimit();

  ~JointPositionLowerLimit();

  // Use default copy constructor.
  JointPositionLowerLimit(const JointPositionLowerLimit&) = default;

  // Use default copy coperator.
  JointPositionLowerLimit& operator=(const JointPositionLowerLimit&) = default;

  // Use default move constructor.
  JointPositionLowerLimit(JointPositionLowerLimit&&) noexcept = default;

  // Use default move assign coperator.
  JointPositionLowerLimit& operator=(JointPositionLowerLimit&&) noexcept 
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
  Eigen::VectorXd qmin_;

};

} // namespace idocp

#endif // IDOCP_JOINT_POSITION_LOWER_LIMIT_HPP_