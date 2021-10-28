#ifndef ROBOTOC_CONSTRAINTS_HPP_
#define ROBOTOC_CONSTRAINTS_HPP_

#include <vector>
#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_solution.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/impulse_constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/constraints/constraints_data.hpp"


namespace robotoc {

///
/// @class Constraints 
/// @brief Stack of the inequality constraints. Composed by constraint 
/// components that inherits ConstraintComponentBase or 
/// ImpulseConstraintComponentBase.
///
class Constraints {
public:
   using ConstraintComponentBasePtr = std::shared_ptr<ConstraintComponentBase>;
   using ImpulseConstraintComponentBasePtr 
      = std::shared_ptr<ImpulseConstraintComponentBase>;

  ///
  /// @brief Default constructor. 
  ///
  Constraints();

  ///
  /// @brief Destructor. 
  ///
  ~Constraints();

  ///
  /// @brief Default copy constructor. 
  ///
  Constraints(const Constraints&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  Constraints& operator=(const Constraints&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  Constraints(Constraints&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  Constraints& operator=(Constraints&&) noexcept = default;

  ///
  /// @brief Appends a constraint component to the cost function.
  /// @param[in] constraint shared pointer to the constraint component appended 
  /// to the constraints.
  ///
  void push_back(const ConstraintComponentBasePtr& constraint);

  ///
  /// @brief Appends a constraint component to the cost function.
  /// @param[in] constraint shared pointer to the constraint component appended 
  /// to the constraints.
  ///
  void push_back(const ImpulseConstraintComponentBasePtr& constraint);

  ///
  /// @brief Clears constraints by removing all components.
  ///
  void clear();

  ///
  /// @brief Checks if the constraints component requres kinematics of robot 
  /// model.
  /// @return true if the constraints component requres kinematics of 
  /// Robot model. false if not.
  ///
  bool useKinematics() const;

  ///
  /// @brief Creates ConstraintsData according to robot model and constraint 
  /// components. 
  /// @param[in] robot Robot model.
  /// @param[in] time_stage Time stage. If -1, the impulse stage is assumed. 
  /// @return Constraints data.
  ///
  ConstraintsData createConstraintsData(const Robot& robot, 
                                        const int time_stage) const;

  ///
  /// @brief Checks whether the current split solution s is feasible or not. 
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data. 
  /// @param[in] s Split solution.
  /// @return true if s is feasible. false if not.
  ///
  bool isFeasible(Robot& robot, ConstraintsData& data,
                  const SplitSolution& s) const;

  ///
  /// @brief Checks whether the current impulse split solution s is feasible or 
  /// not. 
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data. 
  /// @param[in] s Split solution.
  /// @return true if s is feasible. false if not.
  ///
  bool isFeasible(Robot& robot, ConstraintsData& data,
                  const ImpulseSplitSolution& s) const;

  ///
  /// @brief Sets the slack and dual variables of each constraint components. 
  /// @param[in] robot Robot model.
  /// @param[out] data Constraints data. 
  /// @param[in] s Split solution.
  ///
  void setSlackAndDual(Robot& robot, ConstraintsData& data, 
                       const SplitSolution& s) const;

  ///
  /// @brief Sets the slack and dual variables of each impulse constraint 
  /// components. 
  /// @param[in] robot Robot model.
  /// @param[out] data Constraints data. 
  /// @param[in] s Split solution.
  ///
  void setSlackAndDual(Robot& robot, ConstraintsData& data, 
                       const ImpulseSplitSolution& s) const;

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @brief Computes the primal and dual residuals of the constraints. 
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data.
  /// @param[in] s Split solution.
  ///
  void evalConstraint(Robot& robot, ConstraintsData& data, 
                      const SplitSolution& s) const;

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data.
  /// @param[in] s Impulse split solution.
  ///
  void evalConstraint(Robot& robot, ConstraintsData& data, 
                      const ImpulseSplitSolution& s) const;

  ///
  /// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
  /// the products of the Jacobian of the constraints and Lagrange multipliers.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data. 
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual KKT residual.
  ///
  void linearizeConstraints(Robot& robot, ConstraintsData& data, 
                            const double dt, const SplitSolution& s, 
                            SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
  /// the products of the Jacobian of the constraints and Lagrange multipliers.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data. 
  /// @param[in] s Split solution.
  /// @param[out] kkt_residual KKT residual.
  ///
  void linearizeConstraints(Robot& robot, ConstraintsData& data, 
                            const ImpulseSplitSolution& s,
                            ImpulseSplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Condenses the slack and dual variables. linearizeConstraints() must 
  /// be called before this function.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data. 
  /// @param[in] dt Time step.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix Split KKT matrix. The condensed Hessians are added  
  /// to this data.
  /// @param[out] kkt_residual Split KKT residual. The condensed residual are 
  /// added to this data.
  ///
  void condenseSlackAndDual(Robot& robot, ConstraintsData& data,
                            const double dt, const SplitSolution& s,
                            SplitKKTMatrix& kkt_matrix, 
                            SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Condenses the slack and dual variables. linearizeConstraints() must 
  /// be called before this function.
  /// @param[in] robot Robot model.
  /// @param[in] data Constraints data.
  /// @param[in] s Split solution.
  /// @param[out] kkt_matrix Impulse split KKT matrix. The condensed Hessians   
  /// are added to this data.
  /// @param[out] kkt_residual Impulse split KKT residual. The condensed  
  /// residual are added to this data.
  ///
  void condenseSlackAndDual(Robot& robot, ConstraintsData& data,
                            const ImpulseSplitSolution& s,
                            ImpulseSplitKKTMatrix& kkt_matrix, 
                            ImpulseSplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in, out] data Constraints data. 
  /// @param[in] s Split solution.
  /// @param[in] d Split direction.
  ///
  void expandSlackAndDual(ConstraintsData& data, const SplitSolution& s, 
                          const SplitDirection& d) const;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in, out] data Constraints data. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] dts The direction of the switching time regarding of this time 
  /// stage. 
  /// @param[in] s Split solution.
  /// @param[in] d Split direction.
  ///
  void expandSlackAndDual(ConstraintsData& data, const double dt, 
                          const double dts, const SplitSolution& s, 
                          const SplitDirection& d) const;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in, out] data Constraints data. 
  /// @param[in] s Impulse split solution.
  /// @param[in] d Impulse split direction.
  ///
  void expandSlackAndDual(ConstraintsData& data, const ImpulseSplitSolution& s,
                          const ImpulseSplitDirection& d) const;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the directions of the slack variables.
  /// @param[in] data Constraints data.
  /// @return Maximum step size regarding the slack variables.
  ///
  double maxSlackStepSize(const ConstraintsData& data) const;

  ///
  /// @brief Computes and returns the maximum step size by applying 
  /// fraction-to-boundary-rule to the directions of the dual variables.
  /// @param[in] data Constraints data.
  /// @return Maximum step size regarding the dual variables.
  ///
  double maxDualStepSize(const ConstraintsData& data) const;

  ///
  /// @brief Updates the slack variables according to step_size.
  /// @param[in, out] data Constraints data. 
  /// @param[in] step_size Step size. 
  ///
  static void updateSlack(ConstraintsData& data, const double step_size);

  ///
  /// @brief Updates the dual variables according to step_size.
  /// @param[in, out] data Constraints data. 
  /// @param[in] step_size Step size. 
  ///
  static void updateDual(ConstraintsData& data, const double step_size);

  ///
  /// @brief Sets the barrier parameter for all the constraint components.
  /// @param[in] barrier Barrier parameter. Must be positive. Should be small.
  ///
  void setBarrier(const double barrier);

  ///
  /// @brief Sets the parameter of the fraction-to-boundary-rule for all the 
  /// constraint components.
  /// @param[in] fraction_to_boundary_rule Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995.
  ///
  void setFractionToBoundaryRule(const double fraction_to_boundary_rule);

private:
  std::vector<ConstraintComponentBasePtr> position_level_constraints_, 
                                          velocity_level_constraints_, 
                                          acceleration_level_constraints_;
  std::vector<ImpulseConstraintComponentBasePtr> impulse_level_constraints_;

};

} // namespace robotoc

#include "robotoc/constraints/constraints.hxx"

#endif // ROBOTOC_CONSTRAINTS_HPP_