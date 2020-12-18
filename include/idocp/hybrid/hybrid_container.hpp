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
#include "idocp/hybrid/ocp_discretizer.hpp"

#include <vector>
#include <memory>
#include <cassert>

namespace idocp {

template <typename Type, typename ImpulseType>
class hybrid_container;

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
class hybrid_container {
public:
  ///
  /// @brief Construct the standard data, impulse data, and lift data. 
  /// @param[in] robot Robot model.
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data. Default is zero.
  ///
  hybrid_container(const Robot& robot, const int N, const int N_impulse=0) 
    : data(N+1, Type(robot)), 
      aux(N_impulse, Type(robot)), 
      lift(N_impulse, Type(robot)),
      impulse(N_impulse, ImpulseType(robot)) {
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


class OCP {
public:
  ///
  /// @brief Construct only the standard data. 
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] cost Shared ptr to the cost function.
  /// @param[in] constraints Shared ptr to the constraints.
  /// @param[in] T length of the horzion.
  /// @param[in] N number of the standard data.
  /// @param[in] N_impulse number of the impulse data. Default is zero.
  ///
  OCP(const Robot& robot, const std::shared_ptr<CostFunction>& cost, 
      const std::shared_ptr<Constraints>& constraints,
      const double T, const int N, const int N_impulse=0) 
    : data(N, SplitOCP(robot, cost, constraints)), 
      aux(N_impulse, SplitOCP(robot, cost, constraints)), 
      lift(N_impulse, SplitOCP(robot, cost, constraints)),
      impulse(N_impulse, ImpulseSplitOCP(robot, cost->getImpulseCostFunction(), 
                                         constraints->getImpulseConstraints())),
      terminal(TerminalOCP(robot, cost, constraints)),
      discretizer_(T, N, N_impulse) {
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
      discretizer_() {
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

  void discretize(const ContactSequence& contact_sequence, const double t) {
    discretizer_.discretizeOCP(contact_sequence, t);
  }

  const OCPDiscretizer& discretized() const {
    return discretizer_;
  }

  std::vector<SplitOCP> data, aux, lift;
  std::vector<ImpulseSplitOCP> impulse;
  TerminalOCP terminal;

private:
  OCPDiscretizer discretizer_;

};

} // namespace idocp

#endif // IDOCP_HYBRID_CONTAINER_HPP_