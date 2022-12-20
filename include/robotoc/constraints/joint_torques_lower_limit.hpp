#ifndef ROBOTOC_JOINT_TORQUES_LOWER_LIMIT_HPP_
#define ROBOTOC_JOINT_TORQUES_LOWER_LIMIT_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"


namespace robotoc {

///
/// @class JointTorquesLowerLimit
/// @brief Constraint on the lower limits of the joint torques.
///
class JointTorquesLowerLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  JointTorquesLowerLimit(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  JointTorquesLowerLimit();

  ///
  /// @brief Destructor. 
  ///
  ~JointTorquesLowerLimit();

  ///
  /// @brief Default copy constructor. 
  ///
  JointTorquesLowerLimit(const JointTorquesLowerLimit&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  JointTorquesLowerLimit& operator=(const JointTorquesLowerLimit&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JointTorquesLowerLimit(JointTorquesLowerLimit&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JointTorquesLowerLimit& operator=(
      JointTorquesLowerLimit&&) noexcept = default;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const override {}

  bool isFeasible(Robot& robot, const ContactStatus& contact_status, 
                  const GridInfo& grid_info, const SplitSolution& s,
                  ConstraintComponentData& data) const override;

  void setSlack(Robot& robot, const ContactStatus& contact_status, 
                const GridInfo& grid_info, const SplitSolution& s,
                ConstraintComponentData& data) const override;

  void evalConstraint(Robot& robot, const ContactStatus& contact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintComponentData& data) const override;

  void evalDerivatives(Robot& robot, const ContactStatus& contact_status, 
                       const GridInfo& grid_info, const SplitSolution& s,
                       ConstraintComponentData& data,
                       SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(const ContactStatus& contact_status, 
                            const GridInfo& grid_info,
                            ConstraintComponentData& data, 
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(const ContactStatus& contact_status, 
                          const GridInfo& grid_info, const SplitDirection& d, 
                          ConstraintComponentData& data) const override; 

  int dimc() const override;

private:
  int dimc_;
  Eigen::VectorXd umin_;

};

} // namespace robotoc

#endif // ROBOTOC_JOINT_TORQUES_LOWER_LIMIT_HPP_