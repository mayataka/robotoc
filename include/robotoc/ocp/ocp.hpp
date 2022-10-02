#ifndef ROBOTOC_OCP_HPP_
#define ROBOTOC_OCP_HPP_

#include <vector>
#include <memory>
#include <cassert>

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/terminal_ocp.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/sto/sto_cost_function.hpp"
#include "robotoc/sto/sto_constraints.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/time_discretization.hpp"


namespace robotoc {

///
/// @class OCP
/// @brief The (hybrid) optimal control problem.
///
class OCP {
public:
  ///
  /// @brief Construct the optiaml control problem. The switching time 
  /// optimization (STO) algorithm can be enabled with this constructor.
  /// Moreover, the discretization method of the optimal control problem is 
  /// fixed to DiscretizationMethod::PhaseBased.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] sto_cost Shared ptr to the STO cost function.
  /// @param[in] sto_constraints Shared ptr to the STO constraints.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] T Length of the horzion. Must be positive.
  /// @param[in] N Number of the discretization grids of the horizon except for 
  /// the discrete events. Must be positive.
  ///
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints, 
      const std::shared_ptr<STOCostFunction>& sto_cost, 
      const std::shared_ptr<STOConstraints>& sto_constraints, 
      const std::shared_ptr<ContactSequence>& contact_sequence, 
      const double T, const int N);

  ///
  /// @brief Construct the optiaml control problem. The switching time 
  /// optimization (STO) algorithm is disabled with this constructor.
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] T Length of the horzion. Must be positive.
  /// @param[in] N Number of the discretization grids of the horizon except for 
  /// the discrete events. Must be positive.
  ///
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints,  
      const std::shared_ptr<ContactSequence>& contact_sequence,
      const double T, const int N);

  ///
  /// @brief Default Constructor.
  ///
  OCP();

  ///
  /// @brief Default copy constructor. 
  ///
  OCP(const OCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  OCP& operator=(const OCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  OCP(OCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  OCP& operator=(OCP&&) noexcept = default;

  ///
  /// @return const reference to the Robot model. 
  ///
  const Robot& robot() const; 

  ///
  /// @return const reference to the cost function. 
  ///
  const std::shared_ptr<CostFunction>& cost() const;

  ///
  /// @return const reference to the constraints. 
  ///
  const std::shared_ptr<Constraints>& constraints() const;

  ///
  /// @return const reference to the STO cost function. 
  ///
  const std::shared_ptr<STOCostFunction>& sto_cost() const;

  ///
  /// @return const reference to the STO constraints. 
  ///
  const std::shared_ptr<STOConstraints>& sto_constraints() const;

  ///
  /// @return const reference to the STO constraints. 
  ///
  const std::shared_ptr<ContactSequence>& contact_sequence() const;

  ///
  /// @return Length of the horizon. 
  ///
  double T() const;

  ///
  /// @return Number of the discretization grids of the horizon except for 
  /// the discrete events. 
  ///
  int N() const;

  ///
  /// @return Reserved size of the discrete-event data. 
  ///
  int reservedNumDiscreteEvents() const;

  ///
  /// @return true if the switching time optimization (STO) algorithm is enalbed. 
  /// false if not.
  ///
  bool isSTOEnabled() const;

  ///
  /// @brief Split optimal control problem data for the terminal stage.
  ///
  TerminalOCP terminal;

  ///
  /// @brief Displays the optimal control problem onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const OCP& ocp);

private:
  Robot robot_;
  std::shared_ptr<CostFunction> cost_;
  std::shared_ptr<Constraints> constraints_;
  std::shared_ptr<STOCostFunction> sto_cost_;
  std::shared_ptr<STOConstraints> sto_constraints_;
  std::shared_ptr<ContactSequence> contact_sequence_;
  double T_;
  int N_, reserved_num_discrete_events_;
  bool is_sto_enabled_;

  void reserve();

};

} // namespace robotoc

#include "robotoc/ocp/ocp.hxx"

#endif // ROBOTOC_OCP_HPP_ 