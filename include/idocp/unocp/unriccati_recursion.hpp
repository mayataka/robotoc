#ifndef IDOCP_UNRICCATI_RECURSION_HPP_
#define IDOCP_UNRICCATI_RECURSION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"


namespace idocp {

///
/// @class RiccatiRecursion
/// @brief Riccati recursion.
///
class UnRiccatiRecursion {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] nproc Number of the threads in solving the optimal control 
  /// problem. Must be positive. Default is 1.
  ///
  UnRiccatiRecursion(const Robot& robot, const int N);

  ///
  /// @brief Default constructor. 
  ///
  UnRiccatiRecursion();

  ///
  /// @brief Destructor. 
  ///
  ~UnRiccatiRecursion();
 
  ///
  /// @brief Default copy constructor. 
  ///
  UnRiccatiRecursion(const UnRiccatiRecursion&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnRiccatiRecursion& operator=(const UnRiccatiRecursion&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnRiccatiRecursion(UnRiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnRiccatiRecursion& operator=(UnRiccatiRecursion&&) noexcept = default;

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
      const RiccatiFactorization& riccati_factorization_next, 
      UnKKTMatrix& kkt_matrix, UnKKTResidual& kkt_residual, 
      RiccatiFactorization& riccati_factorization);

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
  int N_, nproc_, dimv_;

};

} // namespace idocp

#endif // IDOCP_UNRICCATI_RECURSION_HPP_ 