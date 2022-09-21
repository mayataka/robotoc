#ifndef ROBOTOC_JOINT_POSITION_UPPER_LIMIT_HPP_
#define ROBOTOC_JOINT_POSITION_UPPER_LIMIT_HPP_

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
/// @class JointPositionUpperLimit
/// @brief Constraint on the upper limits of the joint position.
///
class JointPositionUpperLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  ///
  JointPositionUpperLimit(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  JointPositionUpperLimit();

  ///
  /// @brief Destructor. 
  ///
  ~JointPositionUpperLimit();

  ///
  /// @brief Default copy constructor. 
  ///
  JointPositionUpperLimit(const JointPositionUpperLimit&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  JointPositionUpperLimit& operator=(const JointPositionUpperLimit&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  JointPositionUpperLimit(JointPositionUpperLimit&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  JointPositionUpperLimit& operator=(
      JointPositionUpperLimit&&) noexcept = default;

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
  Eigen::VectorXd qmax_;

};

} // namespace robotoc

#endif // ROBOTOC_JOINT_POSITION_UPPER_LIMIT_HPP_