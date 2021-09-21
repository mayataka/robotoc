#ifndef IDOCP_OCP_HPP_
#define IDOCP_OCP_HPP_

#include <vector>
#include <memory>
#include <cassert>

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/impulse_split_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_ocp_discretization.hpp"


namespace idocp {

///
/// @class OCP
/// @brief The (hybrid) optimal control problem.
///
class OCP {
public:
  ///
  /// @brief Construct the optiaml control problem. 
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T Length of the horzion.
  /// @param[in] N Number of the discretization grids.
  /// @param[in] N_impulse Maximum number of the impulses on the horizon. 
  /// Default is 0.
  ///
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints, const double T, 
      const int N, const int N_impulse=0) 
    : data(N, SplitOCP(robot, cost, constraints)), 
      aux(N_impulse, SplitOCP(robot, cost, constraints)), 
      lift(N_impulse, SplitOCP(robot, cost, constraints)),
      impulse(N_impulse, ImpulseSplitOCP(robot, cost, constraints)),
      terminal(TerminalOCP(robot, cost, constraints)),
      time_discretization_(T, N, N_impulse) {
  }

  ///
  /// @brief Default Constructor.
  ///
  OCP() 
    : data(), 
      aux(),
      lift(),
      impulse(),
      terminal(),
      time_discretization_() {
  }

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

  void discretize(const ContactSequence& contact_sequence, const double t) {
    time_discretization_.discretize(contact_sequence, t);
  }

  const HybridOCPDiscretization& discrete() const {
    return time_discretization_;
  }

  std::vector<SplitOCP> data, aux, lift;
  std::vector<ImpulseSplitOCP> impulse;
  TerminalOCP terminal;

private:
  HybridOCPDiscretization time_discretization_;

};

} // namespace idocp

#endif // IDOCP_OCP_HPP_ 