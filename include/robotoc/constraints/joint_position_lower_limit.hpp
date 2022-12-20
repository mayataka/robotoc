#ifndef ROBOTOC_JOINT_POSITION_LOWER_LIMIT_HPP_
#define ROBOTOC_JOINT_POSITION_LOWER_LIMIT_HPP_

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
/// @class JointPositionLowerLimit
/// @brief Constraint on the lower limits of the joint position.
///
class JointPositionLowerLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  JointPositionLowerLimit(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  JointPositionLowerLimit();

  ///
  /// @brief Destructor. 
  ///
  ~JointPositionLowerLimit();

  ///
  /// @brief Default copy constructor. 
  ///
  JointPositionLowerLimit(const JointPositionLowerLimit&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  JointPositionLowerLimit& operator=(const JointPositionLowerLimit&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JointPositionLowerLimit(JointPositionLowerLimit&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JointPositionLowerLimit& operator=(
      JointPositionLowerLimit&&) noexcept = default;

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
  int dimc_, dim_passive_;
  Eigen::VectorXd qmin_;

};

} // namespace robotoc

#endif // ROBOTOC_JOINT_POSITION_LOWER_LIMIT_HPP_