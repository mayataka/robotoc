#ifndef ROBOTOC_CONSTRAINT_COMPONENT_BASE_HPP_
#define ROBOTOC_CONSTRAINT_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @enum KinematicsLevel
/// @brief Kinematics level of the constraint component used in 
/// ConstraintComponentBase.
///
enum class KinematicsLevel {
  PositionLevel,
  VelocityLevel,
  AccelerationLevel
};

///
/// @class ConstraintComponentBase
/// @brief Base class for constraint components. 
///
class ConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995.
  ///
  ConstraintComponentBase(const double barrier, 
                          const double fraction_to_boundary_rule);

  ///
  /// @brief Default constructor. 
  ///
  ConstraintComponentBase();

  ///
  /// @brief Destructor. 
  ///
  virtual ~ConstraintComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  ConstraintComponentBase(const ConstraintComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ConstraintComponentBase& operator=(const ConstraintComponentBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ConstraintComponentBase(ConstraintComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ConstraintComponentBase& operator=(ConstraintComponentBase&&) noexcept 
      = default;

  ///
  /// @brief Checks if the constraint component requres kinematics of robot 
  /// model.
  /// @return true if the constraint component requres kinematics of 
  /// Robot model. false if not.
  ///
  virtual bool useKinematics() const = 0;

  ///
  /// @brief Checks the kinematics level of the constraint component.
  /// @return Kinematics level of the constraint component.
  ///
  virtual KinematicsLevel kinematicsLevel() const = 0;

  ///
  /// @brief Allocates extra data in ConstraintComponentData.
  /// @param[in] data Constraint component data.
  ///
  virtual void allocateExtraData(ConstraintComponentData& data) const = 0;

  ///
  /// @brief Checks whether the current solution s is feasible or not. 
  /// @param[in] robot Robot model.
  /// @param[in] data Constraint data.
  /// @param[in] s Split solution.
  /// @return true if s is feasible. false if not.
  ///
  virtual bool isFeasible(Robot& robot, ConstraintComponentData& data, 
                          const SplitSolution& s) const = 0;

  ///
  /// @brief Sets the slack variables of each constraint components. 
  /// @param[in] robot Robot model.
  /// @param[out] data Constraint data. 
  /// @param[in] s Split solution.
  ///
  virtual void setSlack(Robot& robot, ConstraintComponentData& data, 
                        const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraint data.
  /// @param[in] s Split solution.
  ///
  virtual void evalConstraint(Robot& robot, ConstraintComponentData& data, 
                              const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the derivatives of the priaml residual, i.e., the 
  /// Jacobian of the inequality constraint, and add the product of the 
  /// Jacobian and the dual variable to the KKT residual. This function is 
  /// always called just after evalConstraint().
  /// @param[in] robot Robot model.
  /// @param[in] data Constraint data.
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual Split KKT residual.
  ///
  virtual void evalDerivatives(Robot& robot, ConstraintComponentData& data, 
                               const double dt, const SplitSolution& s, 
                               SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Condenses the slack and dual variables, i.e., factorizes the  
  /// condensed Hessians and KKT residuals. This function is always called 
  /// just after evalDerivatives().
  /// @param[in] data Constraint data.
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix Split KKT matrix. The condensed Hessians are added  
  /// to this object.
  /// @param[out] kkt_residual Split KKT residual. The condensed residuals are 
  /// added to this object.
  ///
  virtual void condenseSlackAndDual(ConstraintComponentData& data,
                                    const double dt, const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix,
                                    SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in, out] data Constraint data.
  /// @param[in] s Split solution.
  /// @param[in] d Split direction.
  ///
  virtual void expandSlackAndDual(ConstraintComponentData& data, 
                                  const SplitSolution& s, 
                                  const SplitDirection& d) const = 0;

  ///
  /// @brief Returns the size of the constraint. 
  /// @return Size of the constraint. 
  /// 
  virtual int dimc() const = 0;

  ///
  /// @brief Sets the slack and dual variables positive.
  /// @param[in, out] data Constraint data.
  ///
  void setSlackAndDualPositive(ConstraintComponentData& data) const;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the direction of the slack variable.
  /// @param[in] data Constraint data.
  /// @return Maximum step size regarding the slack variable.
  ///
  double maxSlackStepSize(const ConstraintComponentData& data) const;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the direction of the dual variable.
  /// @param[in] data Constraint data. 
  /// @return Maximum step size regarding the dual variable.
  ///
  double maxDualStepSize(const ConstraintComponentData& data) const;

  ///
  /// @brief Updates the slack variable according to the step size.
  /// @param[in, out] data Constraint data. 
  /// @param[in] step_size Step size. 
  ///
  static void updateSlack(ConstraintComponentData& data, const double step_size);

  ///
  /// @brief Updates the dual variable according to the step size.
  /// @param[in, out] data Constraint data.
  /// @param[in] step_size Step size. 
  ///
  static void updateDual(ConstraintComponentData& data, const double step_size);

  ///
  /// @brief Returns the barrier parameter.
  ///
  virtual double barrierParameter() const final;

  ///
  /// @brief Returns the parameter of the fraction-to-boundary-rule.
  ///
  virtual double fractionToBoundaryRule() const final;

  ///
  /// @brief Sets the barrier parameter.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  ///
  virtual void setBarrier(const double barrier) final;

  ///
  /// @brief Sets the parameter of the fraction-to-boundary-rule.
  /// @param[in] fraction_to_boundary_rule Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995.
  ///
  virtual void setFractionToBoundaryRule(
      const double fraction_to_boundary_rule) final;

protected:
  ///
  /// @brief Computes the residual in the complementarity slackness between  
  /// the slack and dual variables.
  /// @param[in, out] data Constraint data.
  ///
  void computeComplementarySlackness(ConstraintComponentData& data) const;

  ///
  /// @brief Computes the residual in the complementarity slackness between  
  /// the slack and dual variables.
  /// @param[in, out] data Constraint data.
  /// @param[in] start Start position of the segment.
  /// @param[in] size Size of the segment.
  ///
  void computeComplementarySlackness(ConstraintComponentData& data,
                                     const int start, 
                                     const int size) const;

  ///
  /// @brief Computes the residual in the complementarity slackness between  
  /// the slack and dual variables.
  /// @param[in, out] data Constraint data.
  /// @param[in] start Start position of the segment.
  /// @tparam Size Size of the segment.
  ///
  template <int Size>
  void computeComplementarySlackness(ConstraintComponentData& data,
                                     const int start) const;

  ///
  /// @brief Computes the residual in the complementarity slackness between  
  /// the slack and dual variables.
  /// @param[in] slack An element of the slack variable.
  /// @param[in] dual An element of the dual variable.
  /// @return The complementarity slackness between the slack and dual variables.
  ///
  double computeComplementarySlackness(const double slack, 
                                       const double dual) const;

  ///
  /// @brief Computes the coefficient of the condensing.
  /// @param[in, out] data Constraint component data.
  ///
  static void computeCondensingCoeffcient(ConstraintComponentData& data);

  ///
  /// @brief Computes the coefficient of the condensing.
  /// @param[in, out] data Constraint data.
  /// @param[in] start Start position of the segment.
  /// @param[in] size Size of the segment.
  ///
  static void computeCondensingCoeffcient(ConstraintComponentData& data,
                                          const int start, const int size);

  ///
  /// @brief Computes the coefficient of the condensing.
  /// @param[in, out] data Constraint data.
  /// @param[in] start Start position of the segment.
  /// @tparam Size Size of the segment.
  ///
  template <int Size>
  static void computeCondensingCoeffcient(ConstraintComponentData& data,
                                          const int start);

  ///
  /// @brief Computes the residual in the complementarity slackness between  
  /// the slack and dual variables.
  /// @param[in] slack An element of the slack variable.
  /// @param[in] dual An element of the dual variable.
  /// @param[in] residual An element of the primal residual.
  /// @param[in] cmpl An element of the complementarity slackness.
  /// @return Coefficient of the condensing. 
  ///
  static double computeCondensingCoeffcient(const double slack, 
                                            const double dual, 
                                            const double residual, 
                                            const double cmpl);

  ///
  /// @brief Computes the direction of the dual variable from slack, primal 
  /// residual, complementary slackness, and the direction of the slack.
  /// @param[in, out] data Constraint data.
  ///
  static void computeDualDirection(ConstraintComponentData& data);

  ///
  /// @brief Computes the direction of the dual variable from slack, primal 
  /// residual, complementary slackness, and the direction of the slack.
  /// @param[in, out] data Constraint data.
  /// @param[in] start Start position of the segment.
  /// @param[in] size Size of the segment.
  ///
  static void computeDualDirection(ConstraintComponentData& data, 
                                   const int start, const int size);

  ///
  /// @brief Computes the direction of the dual variable from slack, primal 
  /// residual, complementary slackness, and the direction of the slack.
  /// @param[in, out] data Constraint data.
  /// @param[in] start Start position of the segment.
  /// @tparam Size Size of the segment.
  ///
  template <int Size>
  static void computeDualDirection(ConstraintComponentData& data, 
                                   const int start);

  ///
  /// @brief Computes the direction of the dual variable from slack, primal 
  /// residual, complementary slackness, and the direction of the slack.
  /// @param[in] slack The slack variable.
  /// @param[in] dual The dual variable.
  /// @param[in] dslack The direction of the slack variable.
  /// @param[in] cmpl The residual in the complementary slackness.
  /// @return The direction of the dual variable.
  ///
  static double computeDualDirection(const double slack, const double dual,
                                     const double dslack, const double cmpl);

  ///
  /// @brief Computes the log barrier function of the slack variable.
  /// @param[in] slack Slack variable. All the components must be positive.
  /// @return log barrier function of the slack variable.
  ///
  template <typename VectorType>
  double logBarrier(const Eigen::MatrixBase<VectorType>& slack) const;

private:
  double barrier_, fraction_to_boundary_rule_;

};

} // namespace robotoc

#include "robotoc/constraints/constraint_component_base.hxx"

#endif // ROBOTOC_CONSTRAINT_COMPONENT_BASE_HPP_