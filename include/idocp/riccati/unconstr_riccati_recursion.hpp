#ifndef IDOCP_UNCONSTR_RICCATI_RECURSION_HPP_
#define IDOCP_UNCONSTR_RICCATI_RECURSION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/riccati/unconstr_riccati_factorizer.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"

namespace idocp {

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
  /// @brief Gets of the state feedback gain of the LQR subproblem of the 
  /// specified time stage. 
  /// @param[in] time_stage Time stage of interested. 
  /// @param[in, out] Kq The state feedback gain with respect to the configuration. 
  /// @param[in, out] Kv The state feedback gain with respect to the velocity. 
  ///
  void getStateFeedbackGain(const int time_stage, Eigen::MatrixXd& Kq, 
                            Eigen::MatrixXd& Kv) const;

private:
  int N_;
  double T_, dt_;
  UnconstrRiccatiFactorizer factorizer_;
  std::vector<LQRPolicy> lqr_policy_;

};

} // namespace idocp

#endif // IDOCP_UNCONSTR_RICCATI_RECURSION_HPP_ 