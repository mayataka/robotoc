#ifndef ROBOTOC_RICCATI_RECURSION_HPP_
#define ROBOTOC_RICCATI_RECURSION_HPP_

#include "Eigen/Core"

#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/hybrid/hybrid_container.hpp"
#include "robotoc/riccati/riccati_factorization.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/riccati_factorizer.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"


namespace robotoc {

///
/// @class RiccatiRecursion
/// @brief Riccati recursion solver for hybrid optimal control problems.
/// Solves the KKT system in linear time complexity w.r.t. the length of the 
/// horizon.
///
class RiccatiRecursion {
public:
  ///
  /// @brief Construct a Riccati recursion solver.
  /// @param[in] ocp Optimal control problem. 
  /// @param[in] nthreads Number of the threads used in solving the optimal 
  /// control problem. Must be positive. 
  /// @param[in] max_dts0 Maximum magnitude of the nominal direction of 
  /// the switching time. Used in a heuristic regularization on the dynamic 
  /// programming recursion. Must be positive. Default is 0.1.
  ///
  RiccatiRecursion(const OCP& ocp, const int nthreads, 
                   const double max_dts0=0.1);

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
  /// @brief Sets the regularization on the STO.
  /// @param[in] max_dts0 Maximum magnitude of the nominal direction of 
  /// the switching time. Used in a heuristic regularization on the dynamic 
  /// programming recursion. Must be positive. 
  ///
  void setRegularization(const double max_dts0);

  ///
  /// @brief Reserve the internal data. 
  /// @param[in] ocp Optimal control problem.
  ///
  void reserve(const OCP& ocp);

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] ocp Optimal control problem.
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in, out] factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(const OCP& ocp, KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual, 
                                RiccatiFactorization& factorization);

  ///
  /// @brief Performs the forward Riccati recursion.
  /// @param[in] ocp Optimal control problem.
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] d Direction. d[0].dx is assumed to be computed. 
  ///
  void forwardRiccatiRecursion(const OCP& ocp, const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual, 
                               Direction& d) const;

  ///
  /// @brief Compute the Newton direction in parallel from the Riccati 
  /// factorization factorized by 
  /// RiccatiRecursion::backwardRiccatiRecursion() and 
  /// RiccatiRecursion::forwardRiccatiRecursion().
  /// @param[in] ocp Optimal control problem.
  /// @param[in] contact_sequence Shared ptr to the contact sequence. 
  /// @param[in] factorization Riccati factorization. 
  /// @param[in, out] d Direction. 
  ///
  void computeDirection(OCP& ocp, 
                        const std::shared_ptr<ContactSequence>& contact_sequence, 
                        const RiccatiFactorization& factorization, Direction& d);

  ///
  /// @brief Returns max primal step size.
  /// @return max primal step size.
  /// 
  double maxPrimalStepSize() const;

  ///
  /// @brief Returns max dual step size.
  /// @return max dual step size.
  /// 
  double maxDualStepSize() const;

  ///
  /// @brief Gets of the LQR policies over the horizon. 
  /// @return const reference to the LQR policies.
  ///
  const hybrid_container<LQRPolicy>& getLQRPolicy() const;

private:
  int nthreads_, N_all_;
  RiccatiFactorizer factorizer_;
  hybrid_container<LQRPolicy> lqr_policy_;
  aligned_vector<STOPolicy> sto_policy_;
  SplitRiccatiFactorization factorization_m_;
  Eigen::VectorXd max_primal_step_sizes_, max_dual_step_sizes_;

};

} // namespace robotoc

#endif // ROBOTkeps_OC_RICCATI_RECURSION_HPP_ 