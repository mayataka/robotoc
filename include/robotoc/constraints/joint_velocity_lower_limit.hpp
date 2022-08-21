#ifndef ROBOTOC_JOINT_VELOCITY_LOWER_LIMIT_HPP_
#define ROBOTOC_JOINT_VELOCITY_LOWER_LIMIT_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class JointVelocityLowerLimit
/// @brief Constraint on the lower limits of the joint velocity.
///
class JointVelocityLowerLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  JointVelocityLowerLimit(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  JointVelocityLowerLimit();

  ///
  /// @brief Destructor. 
  ///
  ~JointVelocityLowerLimit();

  ///
  /// @brief Default copy constructor. 
  ///
  JointVelocityLowerLimit(const JointVelocityLowerLimit&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  JointVelocityLowerLimit& operator=(const JointVelocityLowerLimit&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JointVelocityLowerLimit(JointVelocityLowerLimit&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JointVelocityLowerLimit& operator=(
      JointVelocityLowerLimit&&) noexcept = default;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override {}

  bool isFeasible(Robot& robot, const ContactStatus& contact_status, 
                  ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlack(Robot& robot, const ContactStatus& contact_status, 
                ConstraintComponentData& data, 
                const SplitSolution& s) const override;

  void evalConstraint(Robot& robot, const ContactStatus& contact_status, 
                      ConstraintComponentData& data, 
                      const SplitSolution& s) const override;

  void evalDerivatives(Robot& robot, const ContactStatus& contact_status, 
                       ConstraintComponentData& data, const SplitSolution& s,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ContactStatus& contact_status, 
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ContactStatus& contact_status, 
                          ConstraintComponentData& data, 
                          const SplitDirection& d) const override; 

  int dimc() const override;

private:
  int dimc_, dim_passive_;
  Eigen::VectorXd vmin_;

};

} // namespace robotoc

#endif // ROBOTOC_JOINT_VELOCITY_LOWER_LIMIT_HPP_