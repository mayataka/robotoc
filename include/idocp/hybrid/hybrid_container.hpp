#ifndef IDOCP_HYBRID_CONTAINER_HPP_
#define IDOCP_HYBRID_CONTAINER_HPP_

#include "idocp/robot/robot.hpp"

#include "idocp/ocp/split_ocp.hpp"
#include "idocp/impulse/impulse_split_ocp.hpp"
#include "idocp/ocp/terminal_ocp.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/split_riccati_factorizer.hpp"
#include "idocp/impulse/impulse_split_riccati_factorizer.hpp"

#include <vector>
#include <memory>
#include <cassert>

namespace idocp {

template <typename Type, typename ImpulseType>
struct hybrid_container;

using Solution = hybrid_container<SplitSolution, ImpulseSplitSolution>;
using Direction = hybrid_container<SplitDirection, ImpulseSplitDirection>;
using KKTMatrix = hybrid_container<SplitKKTMatrix, ImpulseSplitKKTMatrix>;
using KKTResidual = hybrid_container<SplitKKTResidual, ImpulseSplitKKTResidual>;
using RiccatiFactorization = hybrid_container<SplitRiccatiFactorization, SplitRiccatiFactorization>;
using RiccatiFactorizer = hybrid_container<SplitRiccatiFactorizer, ImpulseSplitRiccatiFactorizer>; 

///
/// @class hybrid_container
/// @brief A container that is useful to formulate the hybrid optimal control 
/// problem. This container has the standard data (with Type), data for lift 
/// stages (with Type), data for aux stages (with Type), and data for impulse 
/// stages (with ImpulseType).
/// @tparam Type The type name of the standard data type.
/// @tparam ImpulseType The type name of the impulse data type.
/// 
///
template <typename Type, typename ImpulseType>
struct hybrid_container {
  ///
  /// @brief Construct the standard data, impulse data, and lift data. 
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  ///
  hybrid_container(const int N, const int N_impulse) 
    : data(N+1, Type()), 
      aux(N_impulse, Type()), 
      lift(N_impulse, Type()),
      impulse(N_impulse, ImpulseType()) {
  }

  ///
  /// @brief Construct the standard data, impulse data, and lift data. 
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  /// @param[in] robot Robot model.
  ///
  hybrid_container(const int N, const int N_impulse, const Robot& robot) 
    : data(N+1, Type(robot)), 
      aux(N_impulse, Type(robot)), 
      lift(N_impulse, Type(robot)),
      impulse(N_impulse, ImpulseType(robot)) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  ///
  hybrid_container(const int N) 
    : data(N+1, Type()), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  /// @param[in] robot Robot model.
  ///
  hybrid_container(const int N, const Robot& robot) 
    : data(N+1, Type(robot)), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Default Constructor.
  ///
  hybrid_container() 
    : data(), 
      aux(),
      lift(),
      impulse() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  hybrid_container(const hybrid_container&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  hybrid_container& operator=(const hybrid_container&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  hybrid_container(hybrid_container&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  hybrid_container& operator=(hybrid_container&&) noexcept = default;

  ///
  /// @brief Overload operator[] to access the standard data, i.e., 
  /// hybrid_container::data as std::vector. 
  ///
  Type& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const Type& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<Type> data, aux, lift;
  std::vector<ImpulseType> impulse;
};


struct OCP {
  ///
  /// @brief Construct the standard data, impulse data, and lift data. 
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  ///
  OCP(const int N, const int N_impulse) 
    : data(N, SplitOCP()), 
      aux(N_impulse, SplitOCP()), 
      lift(N_impulse, SplitOCP()),
      impulse(N_impulse, ImpulseSplitOCP()),
      terminal(TerminalOCP()) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  OCP(const int N, const int N_impulse, const Robot& robot, 
      const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints) 
    : data(N, SplitOCP(robot, cost, constraints)), 
      aux(N_impulse, SplitOCP(robot, cost, constraints)), 
      lift(N_impulse, SplitOCP(robot, cost, constraints)),
      impulse(N_impulse, ImpulseSplitOCP(robot, cost->getImpulseCostFunction(), 
                                         constraints->getImpulseConstraints())),
      terminal(TerminalOCP(robot, cost, constraints)) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  ///
  OCP(const int N) 
    : data(N, SplitOCP()), 
      aux(),
      lift(),
      impulse(),
      terminal(TerminalOCP()) {
  }

  ///
  /// @brief Construct only the standard data. 
  /// @param[in] N number of the standard data.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  ///
  OCP(const int N, const Robot& robot, 
      const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints) 
    : data(N, SplitOCP(robot, cost, constraints)), 
      aux(), 
      lift(),
      impulse(),
      terminal(TerminalOCP(robot, cost, constraints)) {
  }

  ///
  /// @brief Default Constructor.
  ///
  OCP() 
    : data(), 
      aux(),
      lift(),
      impulse(),
      terminal() {
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
  /// hybrid_container::data as std::vector. 
  ///
  SplitOCP& operator[] (const int i) {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  ///
  /// @brief const version of hybrid_container::operator[]. 
  ///
  const SplitOCP& operator[] (const int i) const {
    assert(i >= 0);
    assert(i < data.size());
    return data[i];
  }

  std::vector<SplitOCP> data, aux, lift;
  std::vector<ImpulseSplitOCP> impulse;
  TerminalOCP terminal;
};

} // namespace idocp

#endif // IDOCP_HYBRID_CONTAINER_HPP_