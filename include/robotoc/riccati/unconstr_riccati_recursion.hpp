#ifndef ROBOTOC_UNCONSTR_RICCATI_RECURSION_HPP_
#define ROBOTOC_UNCONSTR_RICCATI_RECURSION_HPP_

#include <vector>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/unconstr/unconstr_ocp.hpp"
#include "robotoc/riccati/unconstr_riccati_factorizer.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"


namespace robotoc {

///
/// @typedef UnconstrRiccatiFactorization
/// @brief Riccati factorization matices of the LQR subproblem for the 
/// unconstrained optimal control problem. 
///
using UnconstrRiccatiFactorization = std::vector<SplitRiccatiFactorization>;

///
/// @class UnconstrRiccatiRecursion
/// @brief Riccati recursion solver for optimal control problems of 
/// unconstrained rigid-body systems.
///
class UnconstrRiccatiRecursion {
public:
  ///
  /// @brief Construct a Riccati recursion solver.
  /// @param[in] ocp Optimial control problem. 
  ///
  UnconstrRiccatiRecursion(const UnconstrOCP& ocp);

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
  /// @param[in, out] factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(KKTMatrix& kkt_matrix, KKTResidual& kkt_residual,
                                UnconstrRiccatiFactorization& factorization);

  ///
  /// @brief Performs the forward Riccati recursion and computes the direction.
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] d Direction. 
  ///
  void forwardRiccatiRecursion(const KKTResidual& kkt_residual, 
                               Direction& d) const;

  ///
  /// @brief Gets of the LQR policies over the horizon. 
  /// @return const reference to the LQR policies.
  ///
  const std::vector<LQRPolicy>& getLQRPolicy() const;

private:
  int N_;
  double dt_;
  UnconstrRiccatiFactorizer factorizer_;
  std::vector<LQRPolicy> lqr_policy_;

};

} // namespace robotoc

#endif // ROBOTOC_UNCONSTR_RICCATI_RECURSION_HPP_ 