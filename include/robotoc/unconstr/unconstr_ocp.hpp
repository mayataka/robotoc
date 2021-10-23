#ifndef ROBOTOC_UNCONSTR_OCP_HPP_
#define ROBOTOC_UNCONSTR_OCP_HPP_

#include <vector>
#include <memory>

#include "robotoc/robot/robot.hpp"

#include "robotoc/unconstr/split_unconstr_ocp.hpp"
#include "robotoc/ocp/terminal_ocp.hpp"
#include "robotoc/cost/cost_function.hpp"
#include "robotoc/constraints/constraints.hpp"


namespace robotoc {

///
/// @class UnconstrOCP
/// @brief An optimal control problem of unconstrained rigid-body systems for 
/// Riccati recursion algorithm.
///
class UnconstrOCP {
public:
  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. 
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] N number of the discretization grids of the horizon.
  ///
  UnconstrOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
              const std::shared_ptr<Constraints>& constraints, const int N) 
    : data(N, SplitUnconstrOCP(robot, cost, constraints)), 
      terminal(TerminalOCP(robot, cost, constraints)) {
  }

  ///
  /// @brief Default Constructor.
  ///
  UnconstrOCP() 
    : data(), 
      terminal() {
  }

  ///
  /// @brief Destructor.
  ///
  ~UnconstrOCP() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrOCP(const UnconstrOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnconstrOCP& operator=(const UnconstrOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrOCP(UnconstrOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrOCP& operator=(UnconstrOCP&&) noexcept = default;

  ///
  /// @brief Overload operator[] to access the data as std::vector. 
  ///
  SplitUnconstrOCP& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const SplitUnconstrOCP& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<SplitUnconstrOCP> data;
  TerminalOCP terminal;
};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_OCP_HPP_ 