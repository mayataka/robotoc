#ifndef ROBOTOC_OCP_HPP_
#define ROBOTOC_OCP_HPP_

#include <vector>
#include <memory>
#include <cassert>

#include "robotoc/robot/robot.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/sto/sto_cost_function.hpp"
#include "robotoc/sto/sto_constraints.hpp"
#include "robotoc/planner/contact_sequence.hpp"


namespace robotoc {

///
/// @class OCP
/// @brief The optimal control problem.
///
struct OCP {
public:
  ///
  /// @brief Construct the optiaml control problem. 
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] sto_cost Shared ptr to the STO cost function.
  /// @param[in] sto_constraints Shared ptr to the STO constraints.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] T Length of the horzion. Must be positive.
  /// @param[in] N Number of the discretization grids of the horizon.
  /// Must be positive.
  /// @param[in] reserved_num_discrete_events Reserved size of the 
  /// discrete-event data. Must be non-negative. Default is 0.
  ///
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints, 
      const std::shared_ptr<STOCostFunction>& sto_cost, 
      const std::shared_ptr<STOConstraints>& sto_constraints, 
      const std::shared_ptr<ContactSequence>& contact_sequence, 
      const double T, const int N, const int reserved_num_discrete_events=0);

  ///
  /// @brief Construct the optiaml control problem. 
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] T Length of the horzion. Must be positive.
  /// @param[in] N Number of the discretization grids of the horizon.
  /// Must be positive.
  /// @param[in] reserved_num_discrete_events Reserved size of the 
  /// discrete-event data. Must be non-negative. Default is 0.
  ///
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints,  
      const std::shared_ptr<ContactSequence>& contact_sequence,
      const double T, const int N, const int reserved_num_discrete_events=0);

  ///
  /// @brief Default constructor.
  ///
  OCP();

  ///
  /// @brief Default destructor.
  ///
  ~OCP() = default;

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
  /// @brief The robot model. 
  ///
  Robot robot; 

  ///
  /// @brief The cost function. 
  ///
  std::shared_ptr<CostFunction> cost = nullptr;

  ///
  /// @brief The constraints. 
  ///
  std::shared_ptr<Constraints> constraints = nullptr;

  ///
  /// @brief The STO cost function. 
  ///
  std::shared_ptr<STOCostFunction> sto_cost= nullptr;

  ///
  /// @brief The STO constraints. 
  ///
  std::shared_ptr<STOConstraints> sto_constraints = nullptr;

  ///
  /// @return The contact sequence. 
  ///
  std::shared_ptr<ContactSequence> contact_sequence = nullptr;

  ///
  /// @return The length of the horizon. 
  ///
  double T = 0.0;

  ///
  /// @return Number of the discretization grids of the horizon.
  ///
  int N = 0;

  ///
  /// @return Reserved size of the discrete-event data. 
  ///
  int reserved_num_discrete_events = 0;

  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const OCP& ocp);
};

} // namespace robotoc

#endif // ROBOTOC_OCP_HPP_ 