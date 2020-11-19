#ifndef IDOCP_RICCATI_RECURSION_HPP_
#define IDOCP_RICCATI_RECURSION_HPP_

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_solution.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/ocp/riccati_factorizer.hpp"
#include "idocp/impulse/impulse_riccati_factorizer.hpp"
#include "idocp/hybrid/hybrid_container.hpp"

namespace idocp {

///
/// @class RiccatiRecursion
/// @brief Riccati factorizer for SplitOCP.
///
class RiccatiRecursion {
public:
  using HybridKKTMatrix = hybrid_container<KKTMatrix, ImpulseKKTMatrix>;
  using HybridKKTResidual = hybrid_container<KKTResidual, ImpulseKKTResidual>;
  using HybridRiccatiFactorization = hybrid_container<RiccatiFactorization, RiccatiFactorization>;
  using HybridRiccatiFactorizer = hybrid_container<RiccatiFactorizer, ImpulseRiccatiFactorizer>;

  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  RiccatiRecursion(const Robot& robot, const double T, const int N, 
                    const int max_num_impulse, const int nproc);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiRecursion();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiRecursion();
 
  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiRecursion(const RiccatiRecursion&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  RiccatiRecursion& operator=(const RiccatiRecursion&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiRecursion(RiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiRecursion& operator=(RiccatiRecursion&&) noexcept = default;

  void backwardRiccatiRecursionTerminal(
      const HybridKKTMatrix& kkt_matrix, const HybridKKTResidual& kkt_residual,
      HybridRiccatiFactorization& riccati_factorization) const;

  void backwardRiccatiRecursion(
      const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
      HybridKKTResidual& kkt_residual, 
      HybridRiccatiFactorization& riccati_factorization);

  void forwardRiccatiRecursionParallel(
      const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
      HybridKKTResidual& kkt_residual);

  void forwardRiccatiRecursionSerial(
      const ContactSequence& contact_sequence, 
      const HybridKKTMatrix& kkt_matrix, const HybridKKTResidual& kkt_residual, 
      HybridRiccatiFactorization& riccati_factorization);

  void backwardStateConstraintFactorization(
      const ContactSequence& contact_sequence, 
      const HybridKKTMatrix& kkt_matrix, 
      std::vector<StateConstraintRiccatiFactorization>& constraint_factorization) const;

private:
  int N_, max_num_impulse_, nproc_;
  double dtau_;
  HybridRiccatiFactorizer riccati_factorizer_;

  void backwardStateConstraintFactorization(
      const ContactSequence& contact_sequence, 
      const HybridKKTMatrix& kkt_matrix, 
      StateConstraintRiccatiFactorization& constraint_factorization,
      const int constraint_index) const;

};

} // namespace idocp

#include "idocp/ocp/riccati_recursion.hxx"

#endif // IDOCP_RICCATI_RECURSION_HPP_ 