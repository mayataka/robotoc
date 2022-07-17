#ifndef ROBOTOC_JOINT_ACCELERATION_LOWER_LIMIT_HPP_
#define ROBOTOC_JOINT_ACCELERATION_LOWER_LIMIT_HPP_

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
/// @class JointAccelerationLowerLimit
/// @brief Constraint on the lower limits of the joint acceleration.
///
class JointAccelerationLowerLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] amin Lower limits of the joint acceleration.
  ///
  JointAccelerationLowerLimit(const Robot& robot, const Eigen::VectorXd& amin);

  ///
  /// @brief Default constructor. 
  ///
  JointAccelerationLowerLimit();

  ///
  /// @brief Destructor. 
  ///
  ~JointAccelerationLowerLimit();

  ///
  /// @brief Default copy constructor. 
  ///
  JointAccelerationLowerLimit(const JointAccelerationLowerLimit&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  JointAccelerationLowerLimit& operator=(
      const JointAccelerationLowerLimit&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JointAccelerationLowerLimit(JointAccelerationLowerLimit&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JointAccelerationLowerLimit& operator=(
      JointAccelerationLowerLimit&&) noexcept = default;

  DEFINE_DEFAULT_CLONE_CONSTRAINT_COMPONENT(JointAccelerationLowerLimit)

  bool useKinematics() const override;

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
  Eigen::VectorXd amin_;

};

} // namespace robotoc

#endif // ROBOTOC_JOINT_ACCELERATION_LOWER_LIMIT_HPP_