#ifndef ROBOTOC_JOINT_VELOCITY_UPPER_LIMIT_HPP_
#define ROBOTOC_JOINT_VELOCITY_UPPER_LIMIT_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class JointVelocityUpperLimit
/// @brief Constraint on the upper limits of the joint velocity.
///
class JointVelocityUpperLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-04.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  JointVelocityUpperLimit(const Robot& robot, const double barrier=1.0e-04,
                          const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Default constructor. 
  ///
  JointVelocityUpperLimit();

  ///
  /// @brief Destructor. 
  ///
  ~JointVelocityUpperLimit();

  ///
  /// @brief Default copy constructor. 
  ///
  JointVelocityUpperLimit(const JointVelocityUpperLimit&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  JointVelocityUpperLimit& operator=(const JointVelocityUpperLimit&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JointVelocityUpperLimit(JointVelocityUpperLimit&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JointVelocityUpperLimit& operator=(
      JointVelocityUpperLimit&&) noexcept = default;

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const {}

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlack(Robot& robot, ConstraintComponentData& data, 
                const SplitSolution& s) const override;

  void evalConstraint(Robot& robot, ConstraintComponentData& data, 
                      const SplitSolution& s) const override;

  void evalDerivatives(Robot& robot, ConstraintComponentData& data, 
                       const double dt, const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const double dt, const SplitSolution& s,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(ConstraintComponentData& data, const SplitSolution& s,
                          const SplitDirection& d) const override; 

  int dimc() const override;

private:
  int dimc_, dim_passive_;
  Eigen::VectorXd vmax_;

};

} // namespace robotoc

#endif // ROBOTOC_JOINT_VELOCITY_UPPER_LIMIT_HPP_