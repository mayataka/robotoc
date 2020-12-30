#ifndef IDOCP_UNRICCATI_RECURSION_HPP_
#define IDOCP_UNRICCATI_RECURSION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/unocp/split_unriccati_factorizer.hpp"

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
  UnRiccatiRecursion(const Robot& robot, const double T, const int N);

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
  /// @param[in] terminal_kkt_matrix KKT matrix at the terminal stage. 
  /// @param[in] terminal_kkt_residual KKT residual at the terminal stage. 
  /// @param[out] riccati_factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursionTerminal(
      const SplitKKTMatrix& terminal_kkt_matrix, 
      const SplitKKTResidual& terminal_kkt_residual,
      RiccatiFactorization& riccati_factorization) const;

  ///
  /// @brief Performs the backward Riccati recursion. Call 
  /// RiccatiRecursion::backwardRiccatiRecursionTerminal() before calling this
  /// function.
  /// @param[in, out] unkkt_matrix KKT matrix. 
  /// @param[in, out] unkkt_residual KKT residual. 
  ///
  void backwardRiccatiRecursion(UnKKTMatrix& unkkt_matrix, 
                                UnKKTResidual& unkkt_residual,
                                RiccatiFactorization& ricccati_factorization);

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] unkkt_residual KKT residual. 
  /// @param[in, out] d Split direction. d[0].dx() must be computed before 
  /// calling this function.
  ///
  void forwardRiccatiRecursion(const UnKKTResidual& unkkt_residual, 
                               Direction& d) const;

private:
  int N_;
  double T_, dtau_;
  std::vector<SplitUnRiccatiFactorizer> factorizer_;

};

} // namespace idocp

#endif // IDOCP_UNRICCATI_RECURSION_HPP_ 