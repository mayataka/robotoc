#ifndef IDOCP_UNCONSTR_OCP_HPP_
#define IDOCP_UNCONSTR_OCP_HPP_

#include <vector>
#include <memory>

#include "idocp/robot/robot.hpp"

#include "idocp/unconstr/split_unconstr_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

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


} // namespace idocp

#endif // IDOCP_UNCONSTR_OCP_HPP_ 