#ifndef IDOCP_UNCONSTR_RICCATI_RECURSION_HPP_
#define IDOCP_UNCONSTR_RICCATI_RECURSION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/riccati/unconstr_riccati_factorizer.hpp"
#include "idocp/riccati/riccati_recursion.hpp"
#include "idocp/riccati/lqr_policy.hpp"

namespace idocp {

///
/// @class UnconstrRiccatiRecursion
/// @brief Riccati recursion solver for optimal control problems of 
/// unconstrained rigid-body systems.
///
class UnconstrRiccatiRecursion {
public:
  ///
  /// @brief Construct a Riccati recursion solver.
  /// @param[in] robot Robot model. 
  /// @param[in] T Length of the horizon. Must be positive.
  /// @param[in] N Number of discretization of the horizon. 
  ///
  UnconstrRiccatiRecursion(const Robot& robot, const double T, const int N);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrRiccatiRecursion();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrRiccatiRecursion();
 
  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrRiccatiRecursion(const UnconstrRiccatiRecursion&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrRiccatiRecursion& operator=(const UnconstrRiccatiRecursion&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrRiccatiRecursion(UnconstrRiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrRiccatiRecursion& operator=(UnconstrRiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[out] riccati_factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(KKTMatrix& unkkt_matrix, 
                                KKTResidual& unkkt_residual,
                                RiccatiFactorization& riccati_factorization);

  ///
  /// @brief Performs the forward Riccati recursion and computes the direction.
  /// @param[in] unkkt_residual KKT residual. 
  /// @param[in, out] d Direction. 
  ///
  void forwardRiccatiRecursion(const KKTResidual& unkkt_residual, 
                               Direction& d) const;

private:
  int N_;
  double T_, dt_;
  UnconstrRiccatiFactorizer factorizer_;
  std::vector<LQRPolicy> lqr_policy_;

};

} // namespace idocp

#endif // IDOCP_UNCONSTR_RICCATI_RECURSION_HPP_ 