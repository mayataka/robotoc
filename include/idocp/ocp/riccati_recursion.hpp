#ifndef IDOCP_RICCATI_RECURSION_HPP_
#define IDOCP_RICCATI_RECURSION_HPP_

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
/// @brief Riccati recursion.
///
class RiccatiRecursion {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of the impulse on the horizon. 
  /// Must be non-negative. Default is 0.
  /// @param[in] nproc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
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

  ///
  /// @brief Performs the backward Riccati recursion for the terminal stage. 
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[out] riccati_factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursionTerminal(
      const HybridKKTMatrix& kkt_matrix, const HybridKKTResidual& kkt_residual,
      HybridRiccatiFactorization& riccati_factorization) const;

  ///
  /// @brief Performs the backward Riccati recursion. Call 
  /// RiccatiRecursion::backwardRiccatiRecursionTerminal() before calling this
  /// function.
  /// @param[in, out] riccati_factorizer Riccati factorizer. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[out] riccati_factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(
      HybridRiccatiFactorizer& riccati_factorizer,
      const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
      HybridKKTResidual& kkt_residual, 
      HybridRiccatiFactorization& riccati_factorization);

  ///
  /// @brief Performs the parallel parts of the forward Riccati recursion. 
  /// Call RiccatiRecursion::backwardRiccatiRecursion() before calling this
  /// function.
  /// @param[in, out] riccati_factorizer Riccati factorizer. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[out] constraint_factorization Pure-state constraint factorization. 
  ///
  void forwardRiccatiRecursionParallel(
      HybridRiccatiFactorizer& riccati_factorizer,
      const ContactSequence& contact_sequence, HybridKKTMatrix& kkt_matrix, 
      HybridKKTResidual& kkt_residual,
      StateConstraintRiccatiFactorization& constraint_factorization);

  ///
  /// @brief Performs the forward factorization of the pure-state constraint
  /// factorization. Call RiccatiRecursion::forwardRiccatiRecursionParallel() 
  /// before calling this function.
  /// @param[in, out] riccati_factorizer Riccati factorizer. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] riccati_factorization Riccati factorization. 
  ///
  void forwardStateConstraintFactorization(
      HybridRiccatiFactorizer& riccati_factorizer,
      const ContactSequence& contact_sequence, 
      const HybridKKTMatrix& kkt_matrix, const HybridKKTResidual& kkt_residual, 
      HybridRiccatiFactorization& riccati_factorization);

  ///
  /// @brief Performs the backward factorization of the pure-state constraint
  /// factorization. Call 
  /// RiccatiRecursion::forwardStateConstraintFactorization() before calling 
  /// this function. 
  /// @param[in] riccati_factorizer Riccati factorizer. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in, out] constraint_factorization Pure-state constraint 
  /// factorization. 
  ///
  void backwardStateConstraintFactorization(
      const HybridRiccatiFactorizer& riccati_factorizer,
      const ContactSequence& contact_sequence, 
      const HybridKKTMatrix& kkt_matrix, 
      StateConstraintRiccatiFactorization& constraint_factorization) const;

  ///
  /// @brief Aggregates the all of the Lagrange multipliers with respect to 
  /// the pure-state constriants for fowrad Riccati recursion.
  /// @param[in] constraint_factorization Pure-state constraint factorization. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] d Split direction.
  /// @param[in, out] riccati_factorization Riccati factorization. 
  ///
  void aggregateLagrangeMultiplierDirection(
      const StateConstraintRiccatiFactorization& constraint_factorization,
      const ContactSequence& contact_sequence, const HybridDirection& d, 
      HybridRiccatiFactorization& riccati_factorization) const;

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] riccati_factorizer Riccati factorizer. 
  /// @param[in] contact_sequence Contact sequence. 
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in] riccati_factorization Riccati factorization. 
  /// @param[in, out] d Split direction. d[0].dx() must be computed before 
  /// calling this function.
  ///
  void forwardRiccatiRecursion(
      const HybridRiccatiFactorizer& riccati_factorizer,
      const ContactSequence& contact_sequence, 
      const HybridKKTMatrix& kkt_matrix, const HybridKKTResidual& kkt_residual, 
      const HybridRiccatiFactorization& riccati_factorization,
      HybridDirection& d);

private:
  int N_, max_num_impulse_, nproc_, dimv_;
  double dtau_;

};

} // namespace idocp

#endif // IDOCP_RICCATI_RECURSION_HPP_ 