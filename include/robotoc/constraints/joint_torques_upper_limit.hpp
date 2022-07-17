#ifndef ROBOTOC_JOINT_TORQUES_UPPER_LIMIT_HPP_
#define ROBOTOC_JOINT_TORQUES_UPPER_LIMIT_HPP_

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
/// @class JointTorquesUpperLimit
/// @brief Constraint on the upper limits of the joint torques.
///
class JointTorquesUpperLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  JointTorquesUpperLimit(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  JointTorquesUpperLimit();

  ///
  /// @brief Destructor. 
  ///
  ~JointTorquesUpperLimit();

  ///
  /// @brief Default copy constructor. 
  ///
  JointTorquesUpperLimit(const JointTorquesUpperLimit&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  JointTorquesUpperLimit& operator=(const JointTorquesUpperLimit&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JointTorquesUpperLimit(JointTorquesUpperLimit&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JointTorquesUpperLimit& operator=(
      JointTorquesUpperLimit&&) noexcept = default;

  DEFINE_DEFAULT_CLONE_CONSTRAINT_COMPONENT(JointTorquesUpperLimit)

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
  int dimc_;
  Eigen::VectorXd umax_;

};

} // namespace robotoc

#endif // ROBOTOC_JOINT_TORQUES_UPPER_LIMIT_HPP_