#ifndef IDOCP_UNCONSTRAINED_CONTAINER_HPP_
#define IDOCP_UNCONSTRAINED_CONTAINER_HPP_

#include "idocp/robot/robot.hpp"

#include "idocp/unocp/split_unocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/unocp/split_unparnmpc.hpp"
#include "idocp/unocp/terminal_unparnmpc.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/unocp/split_unriccati_factorizer.hpp"
#include "idocp/unocp/split_unbackward_correction.hpp"
#include "idocp/hybrid/hybrid_container.hpp"

#include <vector>
#include <memory>
#include <cassert>


namespace idocp {

using UnSolution = std::vector<SplitSolution>;
using UnDirection = std::vector<SplitDirection>;
using UnKKTMatrix = std::vector<SplitUnKKTMatrix>;
using UnKKTResidual = std::vector<SplitUnKKTResidual>;
using UnRiccatiFactorization = std::vector<SplitRiccatiFactorization>;
using UnRiccatiFactorizer = std::vector<SplitUnRiccatiFactorizer>;
using UnBackwardCorrector = std::vector<SplitUnBackwardCorrection>;


struct UnOCP {
public:
  ///
  /// @brief Construct only the standard data. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] N number of the standard data.
  ///
  UnOCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
        const std::shared_ptr<Constraints>& constraints, const int N) 
    : data(N, SplitUnOCP(robot, cost, constraints)), 
      terminal(TerminalOCP(robot, cost, constraints)) {
  }

  ///
  /// @brief Default Constructor.
  ///
  UnOCP() 
    : data(), 
      terminal() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  UnOCP(const UnOCP&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnOCP& operator=(const UnOCP&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnOCP(UnOCP&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnOCP& operator=(UnOCP&&) noexcept = default;

  ///
  /// @brief Overload operator[] to access the standard data, i.e., 
  /// hybrid_container::data as std::vector. 
  ///
  SplitUnOCP& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const SplitUnOCP& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<SplitUnOCP> data;
  TerminalOCP terminal;
};


struct UnParNMPC {
public:
  ///
  /// @brief Construct only the standard data. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] N number of the standard data.
  ///
  UnParNMPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
            const std::shared_ptr<Constraints>& constraints, const int N) 
    : data(N-1, SplitUnParNMPC(robot, cost, constraints)), 
      terminal(TerminalUnParNMPC(robot, cost, constraints)) {
  }

  ///
  /// @brief Default Constructor.
  ///
  UnParNMPC() 
    : data(), 
      terminal() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  UnParNMPC(const UnParNMPC&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  UnParNMPC& operator=(const UnParNMPC&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnParNMPC(UnParNMPC&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnParNMPC& operator=(UnParNMPC&&) noexcept = default;

  ///
  /// @brief Overload operator[] to access the standard data, i.e., 
  /// hybrid_container::data as std::vector. 
  ///
  SplitUnParNMPC& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const SplitUnParNMPC& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<SplitUnParNMPC> data;
  TerminalUnParNMPC terminal;
};

} // namespace idocp

#endif // IDOCP_UNCONSTRAINED_CONTAINER_HPP_ 