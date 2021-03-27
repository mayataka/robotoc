#ifndef IDOCP_UNRICCATI_RECURSION_HPP_
#define IDOCP_UNRICCATI_RECURSION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/unconstrained_container.hpp"
#include "idocp/unocp/split_unriccati_factorizer.hpp"

namespace idocp {

///
/// @class UnRiccatiRecursion
/// @brief Riccati recursion solver for optimal control problems of 
/// unconstrained rigid-body systems.
///
class UnRiccatiRecursion {
public:
  ///
  /// @brief Construct a Riccati recursion solver.
  /// @param[in] robot Robot model. 
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. 
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
      UnRiccatiFactorization& riccati_factorization) const;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in, out] unkkt_matrix KKT matrix. 
  /// @param[in, out] unkkt_residual KKT residual. 
  /// @param[out] riccati_factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(UnKKTMatrix& unkkt_matrix, 
                                UnKKTResidual& unkkt_residual,
                                UnRiccatiFactorization& riccati_factorization);

  ///
  /// @brief Performs the forward Riccati recursion and computes the direction.
  /// @param[in] unkkt_residual KKT residual. 
  /// @param[in, out] d Direction. 
  ///
  void forwardRiccatiRecursion(const UnKKTResidual& unkkt_residual, 
                               UnDirection& d) const;

private:
  int N_;
  double T_, dt_;
  UnRiccatiFactorizer factorizer_;

};

} // namespace idocp

#endif // IDOCP_UNRICCATI_RECURSION_HPP_ 