#ifndef ROBOTOC_IMPULSE_CONSTRAINT_COMPONENT_BASE_HPP_
#define ROBOTOC_IMPULSE_CONSTRAINT_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"


namespace robotoc {

///
/// @class ImpulseConstraintComponentBase
/// @brief Base class for impulse constraint components. 
///
class ImpulseConstraintComponentBase {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-03.
  /// @param[in] fraction_to_boundary_rule Parameter of the 
  /// fraction-to-boundary-rule Must be larger than 0 and smaller than 1. 
  /// Should be between 0.9 and 0.995. Default is 0.995.
  ///
  ImpulseConstraintComponentBase(const double barrier=1.0e-03, 
                                 const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Destructor. 
  ///
  virtual ~ImpulseConstraintComponentBase() {}

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseConstraintComponentBase(
      const ImpulseConstraintComponentBase&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseConstraintComponentBase& operator=(
      const ImpulseConstraintComponentBase&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseConstraintComponentBase(
      ImpulseConstraintComponentBase&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseConstraintComponentBase& operator=(
      ImpulseConstraintComponentBase&&) noexcept = default;

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
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Constraint data.
  /// @param[in] s Impulse split solution.
  /// @return true if s is feasible. false if not.
  ///
  virtual bool isFeasible(Robot& robot, const ImpulseStatus& impulse_status, 
                          ConstraintComponentData& data, 
                          const SplitSolution& s) const = 0;

  ///
  /// @brief Sets the slack variables of each constraint components. 
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[out] data Constraint data. 
  /// @param[in] s Impulse split solution.
  ///
  virtual void setSlack(Robot& robot, const ImpulseStatus& impulse_status, 
                        ConstraintComponentData& data, 
                        const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Constraints data.
  /// @param[in] s Impulse split solution.
  ///
  virtual void evalConstraint(Robot& robot, const ImpulseStatus& impulse_status, 
                              ConstraintComponentData& data, 
                              const SplitSolution& s) const = 0;

  ///
  /// @brief Computes the derivatives of the priaml residual, i.e., the 
  /// Jacobian of the inequality constraint, and add the product of the 
  /// Jacobian and the dual variable to the KKT residual. This function is 
  /// always called just after evalConstraint().
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Constraint data.
  /// @param[in] s Impulse split solution.
  /// @param[out] kkt_residual Impulse split KKT residual.
  ///
  virtual void evalDerivatives(Robot& robot, const ImpulseStatus& impulse_status, 
                               ConstraintComponentData& data, 
                               const SplitSolution& s,
                               SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Condenses the slack and dual variables, i.e., factorizes the  
  /// condensed Hessians and KKT residuals. This function is always called 
  /// just after evalDerivatives().
  /// @param[in] impulse_status Impulse status.
  /// @param[in] data Constraints data.
  /// @param[out] kkt_matrix Impulse split KKT matrix. The condensed Hessians   
  /// are added to this data.
  /// @param[out] kkt_residual Impulse split KKT residual. The condensed KKT
  /// residual are added to this data.
  ///
  virtual void condenseSlackAndDual(
      const ImpulseStatus& impulse_status, ConstraintComponentData& data, 
      SplitKKTMatrix& kkt_matrix, 
      SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] impulse_status Impulse status.
  /// @param[in, out] data Constraints data.
  /// @param[in] d Impulse split direction.
  ///
  virtual void expandSlackAndDual(const ImpulseStatus& impulse_status, 
                                  ConstraintComponentData& data, 
                                  const SplitDirection& d) const = 0;

  ///
  /// @brief Returns the size of the constraint. 
  /// @return Size of the constraints. 
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
  virtual double getBarrierParam() const final;

  ///
  /// @brief Returns the parameter of the fraction-to-boundary-rule.
  ///
  virtual double getFractionToBoundaryRule() const final;

  ///
  /// @brief Sets the barrier parameter.
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  ///
  virtual void setBarrierParam(const double barrier_param) final;

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
  /// residual, complementarity slackness, and the direction of the slack.
  /// @param[in, out] data Constraint data.
  ///
  static void computeDualDirection(ConstraintComponentData& data);

  ///
  /// @brief Computes the direction of the dual variable from slack, primal 
  /// residual, complementarity slackness, and the direction of the slack.
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
  /// @param[in] cmpl The complementary slackness.
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

#include "robotoc/constraints/impulse_constraint_component_base.hxx"

#endif // ROBOTOC_IMPULSE_CONSTRAINT_COMPONENT_BASE_HPP_