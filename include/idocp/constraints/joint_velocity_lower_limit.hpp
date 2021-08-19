#ifndef IDOCP_JOINT_VELOCITY_LOWER_LIMIT_HPP_
#define IDOCP_JOINT_VELOCITY_LOWER_LIMIT_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_base.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

///
/// @class JointVelocityLowerLimit
/// @brief Constraint on the lower limits of the joint velocity.
///
class JointVelocityLowerLimit final : public ConstraintComponentBase {
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
  JointVelocityLowerLimit(const Robot& robot, const double barrier=1.0e-04,
                          const double fraction_to_boundary_rule=0.995);

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

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const {}

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlack(Robot& robot, ConstraintComponentData& data, 
                const SplitSolution& s) const override;

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const SplitSolution& s) const override;

  void computePrimalResidualDerivatives(Robot& robot, ConstraintComponentData& data, 
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
  Eigen::VectorXd vmin_;

};

} // namespace idocp

#endif // IDOCP_JOINT_VELOCITY_LOWER_LIMIT_HPP_