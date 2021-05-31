#ifndef IDOCP_JOINT_TORQUES_LOWER_LIMIT_HPP_
#define IDOCP_JOINT_TORQUES_LOWER_LIMIT_HPP_

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
/// @class JointTorquesLowerLimit
/// @brief Constraint on the lower limits of the joint torques.
///
class JointTorquesLowerLimit final : public ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-04.
  /// @param[in] fraction_to_boundary_rate Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995. Default is 0.995.
  ///
  JointTorquesLowerLimit(const Robot& robot, const double barrier=1.0e-04,
                          const double fraction_to_boundary_rate=0.995);

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

  bool useKinematics() const override;

  KinematicsLevel kinematicsLevel() const override;

  void allocateExtraData(ConstraintComponentData& data) const {}

  bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                  const SplitSolution& s) const override;

  void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                       const SplitSolution& s) const override;

  void augmentDualResidual(Robot& robot, ConstraintComponentData& data, 
                           const double dt, const SplitSolution& s,
                           SplitKKTResidual& kkt_residual) const override;

  void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                            const double dt, const SplitSolution& s,
                            SplitKKTMatrix& kkt_matrix,
                            SplitKKTResidual& kkt_residual) const override;

  void expandSlackAndDual(ConstraintComponentData& data, const SplitSolution& s,
                          const SplitDirection& d) const override; 

  void computePrimalAndDualResidual(Robot& robot, ConstraintComponentData& data, 
                                    const SplitSolution& s) const override;

  int dimc() const override;

private:
  int dimc_;
  Eigen::VectorXd umin_;

};

} // namespace idocp

#endif // IDOCP_JOINT_TORQUES_LOWER_LIMIT_HPP_