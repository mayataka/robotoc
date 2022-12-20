#ifndef ROBOTOC_CONSTRAINTS_HPP_
#define ROBOTOC_CONSTRAINTS_HPP_

#include <vector>
#include <memory>
#include <unordered_map>
#include <iostream>

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/constraints/constraint_component_base.hpp"
#include "robotoc/constraints/impact_constraint_component_base.hpp"
#include "robotoc/constraints/constraint_component_data.hpp"
#include "robotoc/constraints/constraints_data.hpp"
#include "robotoc/ocp/grid_info.hpp"


namespace robotoc {

///
/// @class Constraints 
/// @brief Stack of the inequality constraints. Composed by constraint 
/// components that inherits ConstraintComponentBase or 
/// ImpactConstraintComponentBase.
///
class Constraints {
public:
   using ConstraintComponentBasePtr = std::shared_ptr<ConstraintComponentBase>;
   using ImpactConstraintComponentBasePtr 
      = std::shared_ptr<ImpactConstraintComponentBase>;

  ///
  /// @brief Constructor. 
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  /// Default is 1.0e-03
  /// @param[in] fraction_to_boundary_rule Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995. Default is 0.995.
  ///
  Constraints(const double barrier_param=1.0e-03, 
              const double fraction_to_boundary_rule=0.995);

  ///
  /// @brief Default destructor. 
  ///
  ~Constraints() = default;

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
  /// @brief Checks if thsi has a constraint component of the specified name. 
  /// @param[in] name Name of the constraint component.
  /// @return treu if a constraint component of the specified name exists. 
  ///
  bool exist(const std::string& name) const;

  ///
  /// @brief Adds a constraint component. If a component of the same name 
  /// exists, throws an exeption.
  /// @param[in] name Name of the constraint component.
  /// @param[in] constraint shared pointer to the constraint component.
  ///
  void add(const std::string& name, ConstraintComponentBasePtr constraint);

  ///
  /// @brief Adds a constraint component. If a component of the same name 
  /// exists, throws an exeption.
  /// @param[in] name Name of the constraint component.
  /// @param[in] constraint shared pointer to the constraint component.
  ///
  void add(const std::string& name, ImpactConstraintComponentBasePtr constraint);

  ///
  /// @brief Erases a constraint component. If a component of the specified 
  /// name does not exist, throws an exeption.
  /// @param[in] name Name of the constraint component.
  ///
  void erase(const std::string& name);

  ///
  /// @brief Gets a constraint component. If a component of the specified 
  /// name does not exist, throws an exeption. 
  /// @param[in] name Name of the constraint component.
  /// @return Shared ptr to the specified constraint component.
  ///
  ConstraintComponentBasePtr getConstraintComponent(const std::string& name) const;

  ///
  /// @brief Gets a constraint component. If a component of the specified 
  /// name does not exist, throws an exeption. 
  /// @param[in] name Name of the constraint component.
  /// @return Shared ptr to the specified constraint component.
  ///
  ImpactConstraintComponentBasePtr getImpactConstraintComponent(const std::string& name) const;

  ///
  /// @brief Clears constraints by removing all components.
  ///
  void clear();

  ///
  /// @brief Creates ConstraintsData according to robot model and constraint 
  /// components. 
  /// @param[in] robot Robot model.
  /// @param[in] time_stage Time stage. If -1, the impact stage is assumed. 
  /// @return Constraints data.
  ///
  ConstraintsData createConstraintsData(const Robot& robot, 
                                        const int time_stage=-1) const;

  ///
  /// @brief Checks whether the current split solution s is feasible or not. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] data Constraints data. 
  /// @return true if s is feasible. false if not.
  ///
  bool isFeasible(Robot& robot, const ContactStatus& contact_status,
                  const GridInfo& grid_info, const SplitSolution& s,
                  ConstraintsData& data) const;

  ///
  /// @brief Checks whether the current impact split solution s is feasible or 
  /// not. 
  /// @param[in] robot Robot model.
  /// @param[in] impact_status Impact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] data Constraints data. 
  /// @return true if s is feasible. false if not.
  ///
  bool isFeasible(Robot& robot, const ImpactStatus& impact_status, 
                  const GridInfo& grid_info, const SplitSolution& s,
                  ConstraintsData& data) const;

  ///
  /// @brief Sets the slack and dual variables of each constraint components. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] data Constraints data. 
  ///
  void setSlackAndDual(Robot& robot, const ContactStatus& contact_status, 
                       const GridInfo& grid_info, const SplitSolution& s,
                       ConstraintsData& data) const;

  ///
  /// @brief Sets the slack and dual variables of each impact constraint 
  /// components. 
  /// @param[in] robot Robot model.
  /// @param[in] impact_status Impact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] data Constraints data. 
  ///
  void setSlackAndDual(Robot& robot, const ImpactStatus& impact_status, 
                       const GridInfo& grid_info, const SplitSolution& s,
                       ConstraintsData& data) const;

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @brief Computes the primal and dual residuals of the constraints. 
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] data Constraints data.
  ///
  void evalConstraint(Robot& robot, const ContactStatus& contact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintsData& data) const;

  ///
  /// @brief Computes the primal residual, residual in the complementary 
  /// slackness, and the log-barrier function of the slack varible.
  /// @param[in] robot Robot model.
  /// @param[in] impact_status Impact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Impact split solution.
  /// @param[in, out] data Constraints data.
  ///
  void evalConstraint(Robot& robot, const ImpactStatus& impact_status, 
                      const GridInfo& grid_info, const SplitSolution& s,
                      ConstraintsData& data) const;

  ///
  /// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
  /// the products of the Jacobian of the constraints and Lagrange multipliers.
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] data Constraints data. 
  /// @param[in, out] kkt_residual KKT residual.
  ///
  void linearizeConstraints(Robot& robot, const ContactStatus& contact_status, 
                            const GridInfo& grid_info, const SplitSolution& s, 
                            ConstraintsData& data,
                            SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Evaluates the constraints (i.e., calls evalConstraint()) and adds 
  /// the products of the Jacobian of the constraints and Lagrange multipliers.
  /// @param[in] robot Robot model.
  /// @param[in] impact_status Impact status.
  /// @param[in] grid_info Grid info.
  /// @param[in] s Split solution.
  /// @param[in, out] data Constraints data. 
  /// @param[in, out] kkt_residual KKT residual.
  ///
  void linearizeConstraints(Robot& robot, const ImpactStatus& impact_status, 
                            const GridInfo& grid_info, const SplitSolution& s, 
                            ConstraintsData& data,
                            SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Condenses the slack and dual variables. linearizeConstraints() must 
  /// be called before this function.
  /// @param[in] contact_status Contact status.
  /// @param[in] grid_info Grid info.
  /// @param[in, out] data Constraints data. 
  /// @param[in, out] kkt_matrix Split KKT matrix. The condensed Hessians are added  
  /// to this data.
  /// @param[in, out] kkt_residual Split KKT residual. The condensed residual are 
  /// added to this data.
  ///
  void condenseSlackAndDual(const ContactStatus& contact_status, 
                            const GridInfo& grid_info, ConstraintsData& data,
                            SplitKKTMatrix& kkt_matrix, 
                            SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Condenses the slack and dual variables. linearizeConstraints() must 
  /// be called before this function.
  /// @param[in] impact_status Impact status.
  /// @param[in] grid_info Grid info.
  /// @param[in, out] data Constraints data.
  /// @param[in, out] kkt_matrix Impact split KKT matrix. The condensed Hessians   
  /// are added to this data.
  /// @param[in, out] kkt_residual Impact split KKT residual. The condensed  
  /// residual are added to this data.
  ///
  void condenseSlackAndDual(const ImpactStatus& impact_status, 
                            const GridInfo& grid_info, ConstraintsData& data, 
                            SplitKKTMatrix& kkt_matrix, 
                            SplitKKTResidual& kkt_residual) const;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] contact_status Contact status.
  /// @param[in] grid_info Grid info.
  /// @param[in, out] data Constraint data.
  /// @param[in] d Split direction.
  ///
  void expandSlackAndDual(const ContactStatus& contact_status,
                          const GridInfo& grid_info, const SplitDirection& d,
                          ConstraintsData& data) const;

  ///
  /// @brief Expands the slack and dual, i.e., computes the directions of the 
  /// slack and dual variables from the directions of the primal variables.
  /// @param[in] impact_status Impact status.
  /// @param[in] grid_info Grid info.
  /// @param[in, out] data Constraint data.
  /// @param[in] d Split direction.
  ///
  void expandSlackAndDual(const ImpactStatus& impact_status,
                          const GridInfo& grid_info, const SplitDirection& d,
                          ConstraintsData& data) const;

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
  /// @param[in] barrier_param Barrier parameter. Must be positive. Should be small.
  ///
  void setBarrierParam(const double barrier_param);

  ///
  /// @brief Sets the parameter of the fraction-to-boundary-rule for all the 
  /// constraint components.
  /// @param[in] fraction_to_boundary_rule Must be larger than 0 and smaller 
  /// than 1. Should be between 0.9 and 0.995.
  ///
  void setFractionToBoundaryRule(const double fraction_to_boundary_rule);

  ///
  /// @brief Gets the barrier parameter.
  /// @return Barrier parameter. 
  ///
  double getBarrierParam() const;

  ///
  /// @brief Gets the parameter of the fraction-to-boundary-rule. 
  /// @return The parameter of the fraction-to-boundary-rule. 
  ///
  double getFractionToBoundaryRule() const;

  ///
  /// @brief Gets a list of the position-level constraints. 
  /// @return Name list of the position-level constraints.
  ///
  std::vector<std::string> getPositionLevelConstraintList() const;

  ///
  /// @brief Gets a list of the velocity-level constraints. 
  /// @return Name list of the velocity-level constraints.
  ///
  std::vector<std::string> getVelocityLevelConstraintList() const;

  ///
  /// @brief Gets a list of the acceleration-level constraints. 
  /// @return Name list of the acceleration-level constraints.
  ///
  std::vector<std::string> getAccelerationLevelConstraintList() const;

  ///
  /// @brief Gets a list of the impact-level constraints. 
  /// @return Name list of acceleration impact-level constraints.
  ///
  std::vector<std::string> getImpactLevelConstraintList() const;

  ///
  /// @brief Gets a list of the constraints. 
  /// @return Name list of acceleration impact-level constraints.
  ///
  std::vector<std::string> getConstraintList() const;

  ///
  /// @brief Displays the constraints onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const Constraints& constraints);

  friend std::ostream& operator<<(std::ostream& os, 
                                  const std::shared_ptr<Constraints>& constraints);

private:
  std::vector<ConstraintComponentBasePtr> position_level_constraints_, 
                                          velocity_level_constraints_, 
                                          acceleration_level_constraints_;
  std::vector<ImpactConstraintComponentBasePtr> impact_level_constraints_;
  std::unordered_map<std::string, size_t> position_level_constraint_names_, 
                                          velocity_level_constraint_names_, 
                                          acceleration_level_constraint_names_,
                                          impact_level_constraint_names_;
  double barrier_, fraction_to_boundary_rule_;
};

} // namespace robotoc

#endif // ROBOTOC_CONSTRAINTS_HPP_