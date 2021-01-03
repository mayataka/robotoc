#ifndef IDOCP_RICCATI_RECURSION_HPP_
#define IDOCP_RICCATI_RECURSION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/ocp/state_constraint_riccati_factorization.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/hybrid/ocp_discretizer.hpp"

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
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nthreads Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  RiccatiRecursion(const Robot& robot, const int N, const int nthreads);

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
      const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
      RiccatiFactorization& riccati_factorization) const;

  ///
  /// @brief Performs the backward Riccati recursion. Call 
  /// RiccatiRecursion::backwardRiccatiRecursionTerminal() before calling this
  /// function.
  /// @param[in, out] riccati_factorizer Riccati factorizer. 
  /// @param[in] ocp_discretizer OCP discretizer.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[out] riccati_factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(
      RiccatiFactorizer& riccati_factorizer, 
      const OCPDiscretizer& ocp_discretizer, KKTMatrix& kkt_matrix, 
      KKTResidual& kkt_residual, RiccatiFactorization& riccati_factorization);

  ///
  /// @brief Performs the parallel parts of the forward Riccati recursion. 
  /// Call RiccatiRecursion::backwardRiccatiRecursion() before calling this
  /// function.
  /// @param[in, out] riccati_factorizer Riccati factorizer. 
  /// @param[in] ocp_discretizer OCP discretizer.
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[out] constraint_factorization Pure-state constraint factorization. 
  ///
  void forwardRiccatiRecursionParallel(
      RiccatiFactorizer& riccati_factorizer, 
      const OCPDiscretizer& ocp_discretizer, KKTMatrix& kkt_matrix, 
      KKTResidual& kkt_residual,
      StateConstraintRiccatiFactorization& constraint_factorization);

  ///
  /// @brief Performs the forward factorization of the pure-state constraint
  /// factorization. Call RiccatiRecursion::forwardRiccatiRecursionParallel() 
  /// before calling this function.
  /// @param[in, out] riccati_factorizer Riccati factorizer. 
  /// @param[in] ocp_discretizer OCP discretizer.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] riccati_factorization Riccati factorization. 
  ///
  void forwardStateConstraintFactorization(
      RiccatiFactorizer& riccati_factorizer, 
      const OCPDiscretizer& ocp_discretizer, const KKTMatrix& kkt_matrix, 
      const KKTResidual& kkt_residual, 
      RiccatiFactorization& riccati_factorization);

  ///
  /// @brief Performs the backward factorization of the pure-state constraint
  /// factorization. Call 
  /// RiccatiRecursion::forwardStateConstraintFactorization() before calling 
  /// this function. 
  /// @param[in] riccati_factorizer Riccati factorizer. 
  /// @param[in] ocp_discretizer OCP discretizer.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in, out] constraint_factorization Pure-state constraint 
  /// factorization. 
  ///
  void backwardStateConstraintFactorization(
      const RiccatiFactorizer& riccati_factorizer, 
      const OCPDiscretizer& ocp_discretizer, const KKTMatrix& kkt_matrix, 
      StateConstraintRiccatiFactorization& constraint_factorization) const;

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] riccati_factorizer Riccati factorizer. 
  /// @param[in] ocp_discretizer OCP discretizer.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in] riccati_factorization Riccati factorization. 
  /// @param[in, out] d Split direction. d[0].dx() must be computed before 
  /// calling this function.
  ///
  void forwardRiccatiRecursion(
      const RiccatiFactorizer& riccati_factorizer,
      const OCPDiscretizer& ocp_discretizer, const KKTMatrix& kkt_matrix, 
      const KKTResidual& kkt_residual, 
      const RiccatiFactorization& riccati_factorization, Direction& d);

private:
  int N_, nthreads_, dimv_;

};

} // namespace idocp

#endif // IDOCP_RICCATI_RECURSION_HPP_ 