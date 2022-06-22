#ifndef ROBOTOC_OCP_HPP_
#define ROBOTOC_OCP_HPP_

#include <vector>
#include <memory>
#include <cassert>

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_ocp.hpp"
#include "robotoc/impulse/impulse_split_ocp.hpp"
#include "robotoc/ocp/terminal_ocp.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"
#include "robotoc/hybrid/sto_cost_function.hpp"
#include "robotoc/hybrid/sto_constraints.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/hybrid/time_discretization.hpp"


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
  /// @brief Resize the internal data without reallocating all the data. 
  /// @param[in] max_num_each_discrete_events Maximum possible number of the 
  /// each discrete events on the horizon. Must be non-negative. Default is 0.
  ///
  void resize(const int max_num_each_discrete_events);

  ///
  /// @brief Sets the discretization method of the optimal contro problem. 
  /// @param[in] discretization_method The discretization method.
  ///
  void setDiscretizationMethod(const DiscretizationMethod discretization_method);

  ///
  /// @brief Discretizes the optimal control problem according to the 
  /// input current contact sequence and intial time of the horizon.
  /// @param[in] t Initial time of the horizon. 
  ///
  void discretize(const double t);

  ///
  /// @brief Performs mesh-refimenent while keeping the total number of the 
  /// discretization grids over the horizon. This function does it only when
  /// the discretization method is set to DiscretizationMethod::PhaseBased.
  /// Otherwise, does nothing.
  /// @param[in] t Initial time of the horizon. 
  ///
  void meshRefinement(const double t);

  ///
  /// @brief Returns the discrete-time formulation, that is, the discretization 
  /// of the optimal control problem. 
  /// @return The discretization of the optimal control problem. 
  ///
  const TimeDiscretization& discrete() const;

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
  /// @return Maximum possible number of the each discrete events on the horizon. 
  ///
  int maxNumEachDiscreteEvents() const;

  ///
  /// @return true if the switching time optimization (STO) algorithm is enalbed. 
  /// false if not.
  ///
  bool isSTOEnabled() const;

  ///
  /// @brief Overload operator[] to access the standard data, i.e., 
  /// OCP::data as std::vector. 
  ///
  SplitOCP& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of OCP::operator[]. 
  ///
  const SplitOCP& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief Split optimal control problem data for the time stages.
  ///
  std::vector<SplitOCP> data;

  ///
  /// @brief Split optimal control problem data for the auxiliary stages 
  /// (additional time stages just after the impulse events).
  ///
  std::vector<SplitOCP> aux;

  ///
  /// @brief Split optimal control problem data for the lift stages 
  /// (additional time stages just after the lift events).
  ///
  std::vector<SplitOCP> lift;

  ///
  /// @brief Split optimal control problem data for the impulse stages 
  /// (additional time stages at the impulse events).
  ///
  std::vector<ImpulseSplitOCP> impulse;

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
  TimeDiscretization discretization_;
  double T_;
  int N_, max_num_each_discrete_events_;
  bool is_sto_enabled_;

};

} // namespace robotoc

#include "robotoc/ocp/ocp.hxx"

#endif // ROBOTOC_OCP_HPP_ 