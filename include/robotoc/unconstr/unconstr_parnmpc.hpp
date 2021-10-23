#ifndef ROBOTOC_UNCONSTR_PARNMPC_HPP_
#define ROBOTOC_UNCONSTR_PARNMPC_HPP_

#include <vector>
#include <memory>

#include "robotoc/robot/robot.hpp"

#include "robotoc/unconstr/split_unconstr_parnmpc.hpp"
#include "robotoc/unconstr/terminal_unconstr_parnmpc.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"


namespace robotoc {

///
/// @class UnconstrParNMPC
/// @brief An optimal control problem of unconstrained rigid-body systems for 
/// ParNMPC algorithm.
///
class UnconstrParNMPC {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] N number of the discretization grids of the horizon.
  ///
  UnconstrParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
                  const std::shared_ptr<Constraints>& constraints, const int N) 
    : data(N-1, SplitUnconstrParNMPC(robot, cost, constraints)), 
      terminal(TerminalUnconstrParNMPC(robot, cost, constraints)) {
  }

  ///
  /// @brief Default Constructor.
  ///
  UnconstrParNMPC() 
    : data(), 
      terminal() {
  }

  ///
  /// @brief Destructor.
  ///
  ~UnconstrParNMPC() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrParNMPC(const UnconstrParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnconstrParNMPC& operator=(const UnconstrParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrParNMPC(UnconstrParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrParNMPC& operator=(UnconstrParNMPC&&) noexcept = default;

  ///
  /// @brief Overload operator[] to access the data as std::vector. 
  ///
  SplitUnconstrParNMPC& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const SplitUnconstrParNMPC& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<SplitUnconstrParNMPC> data;
  TerminalUnconstrParNMPC terminal;
};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_PARNMPC_HPP_ 