#ifndef IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_
#define IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/constraints/constraint_component_data.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"


namespace idocp {

///
/// @class KinematicsLevel
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
  /// @param[in] fraction_to_boundary_rate Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995.
  ///
  ConstraintComponentBase(const double barrier, 
                          const double fraction_to_boundary_rate);

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
  /// @brief Sets the slack and dual variables of each constraint components. 
  /// @param[in] robot Robot model.
  /// @param[out] data Constraint data. 
  /// @param[in] s Split solution.
  ///
  virtual void setSlackAndDual(Robot& robot, ConstraintComponentData& data, 
                               const SplitSolution& s) const = 0;

  ///
  /// @brief Augments the dual residual of the constraint to the KKT residual.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraint data.
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual Split KKT residual.
  ///
  virtual void augmentDualResidual(Robot& robot, ConstraintComponentData& data,
                                   const double dt, const SplitSolution& s,
                                   SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Condenses the slack and dual variables and factorizes the condensed 
  /// Hessians and KKT residuals.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraint data.
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix Split KKT matrix. The condensed Hessians are added  
  /// to this object.
  /// @param[out] kkt_residual Split KKT residual. The condensed residuals are 
  /// added to this object.
  ///
  virtual void condenseSlackAndDual(Robot& robot, ConstraintComponentData& data,
                                    const double dt, const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix,
                                    SplitKKTResidual& kkt_residual) const = 0;

  ///
  /// @brief Computes the directions of the slack and dual variables.
  /// @param[in] robot Robot model.
  /// @param[in, out] data Constraint data.
  /// @param[in] s Split solution.
  /// @param[in] d Split direction.
  ///
  virtual void computeSlackAndDualDirection(
      Robot& robot, ConstraintComponentData& data, const SplitSolution& s, 
      const SplitDirection& d) const = 0;

  ///
  /// @brief Computes the primal and dual residuals of the constraint. 
  /// @param[in] robot Robot model.
  /// @param[in] data Constraint data.
  /// @param[in] s Split solution.
  ///
  virtual void computePrimalAndDualResidual(
      Robot& robot, ConstraintComponentData& data, 
      const SplitSolution& s) const = 0;

  ///
  /// @brief Returns the size of the constraint. 
  /// @return Size of the constraint. 
  /// 
  virtual int dimc() const = 0;

  ///
  /// @brief Returns the l1-norm of the primal residual of the constraint. 
  /// Before calling this function, 
  /// ConstraintComponentBase::computePrimalResidual() must be called.
  /// @param[in] data Constraint data. 
  /// @return l1-norm of the primal residual of the constraint. 
  ///
  virtual double l1NormPrimalResidual(
      const ConstraintComponentData& data) const final;

  ///
  /// @brief Returns the squared norm of the primal and dual residuals of the 
  /// constraint. Before calling this function, 
  /// ConstraintComponentBase::computePrimalResidual() must be called.
  /// @param[in] data Constraint data.
  /// @return Squared norm of the primal and dual residuals of the constraint. 
  ///
  virtual double squaredNormPrimalAndDualResidual(
      const ConstraintComponentData& data) const final;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the direction of the slack variable.
  /// @param[in] data Constraint data.
  /// @return Maximum step size regarding the slack variable.
  ///
  virtual double maxSlackStepSize(
      const ConstraintComponentData& data) const final;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the direction of the dual variable.
  /// @param[in] data Constraint data. 
  /// @return Maximum step size regarding the dual variable.
  ///
  virtual double maxDualStepSize(
      const ConstraintComponentData& data) const final;

  ///
  /// @brief Updates the slack variable according to the step size.
  /// @param[in, out] data Constraint data. 
  /// @param[in] step_size Step size. 
  ///
  virtual void updateSlack(ConstraintComponentData& data, 
                           const double step_size) const final;

  ///
  /// @brief Updates the dual variable according to the step size.
  /// @param[in, out] data Constraint data.
  /// @param[in] step_size Step size. 
  ///
  virtual void updateDual(ConstraintComponentData& data, 
                          const double step_size) const final;

  ///
  /// @brief Computes and returns the value of the barrier function of the slack 
  /// variable.
  /// @param[in] data Constraint data. 
  /// @return Value of the barrier function. 
  ///
  virtual double costSlackBarrier(
      const ConstraintComponentData& data) const final;

  ///
  /// @brief Computes and returns the value of the barrier function of the slack 
  /// variable with the step size.
  /// @param[in] data Constraint data.
  /// @param[in] step_size Step size. 
  /// @return Value of the barrier function. 
  ///
  virtual double costSlackBarrier(const ConstraintComponentData& data, 
                                  const double step_size) const final;

  ///
  /// @brief Returns the barrier parameter.
  ///
  virtual double barrier() const final;

  ///
  /// @brief Returns the rate of the fraction-to-boundary-rule.
  ///
  virtual double fractionToBoundaryRate() const final;

  ///
  /// @brief Sets the barrier parameter.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  ///
  virtual void setBarrier(const double barrier) final;

  ///
  /// @brief Sets the fraction to boundary rate.
  /// @param[in] fraction_to_boundary_rate Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995.
  ///
  virtual void setFractionToBoundaryRate(
      const double fraction_to_boundary_rate) final;

protected:
  ///
  /// @brief Sets the slack and dual variables positive.
  /// @param[in, out] data Constraint data.
  ///
  virtual void setSlackAndDualPositive(
      ConstraintComponentData& data) const final;

  ///
  /// @brief Computes the duality residual between the slack and dual variables.
  /// @param[in, out] data Constraint data.
  ///
  virtual void computeDuality(ConstraintComponentData& data) const final;

  ///
  /// @brief Computes the direction of the dual variable from slack, residual,
  /// duality, and the direction of the slack.
  /// @param[in, out] data Constraint data.
  ///
  virtual void computeDualDirection(ConstraintComponentData& data) const final;

  ///
  /// @brief Computes the duality residual between the slack and dual variables.
  /// @param[in] slack An element of the slack variable.
  /// @param[in] dual An element of the dual variable.
  /// @return An element of the duality of the slack and dual variables.
  ///
  virtual double computeDuality(const double slack, 
                                const double dual) const final;

  ///
  /// @brief Computes the direction of the dual variable from slack, residual,
  /// duality, and the direction of the slack.
  /// @param[in] slack An element of the slack variable.
  /// @param[in] dual An element of the dual variable.
  /// @param[in] dslack An element of the direction of the slack variable.
  /// @param[in] duality An element of the duality.
  /// @return An element of the direction of the dual variable.
  ///
  virtual double computeDualDirection(const double slack, const double dual,
                                      const double dslack, 
                                      const double duality) const final;

private:
  double barrier_, fraction_to_boundary_rate_;

};

} // namespace idocp

#include "idocp/constraints/constraint_component_base.hxx"

#endif // IDOCP_CONSTRAINT_COMPONENT_BASE_HPP_