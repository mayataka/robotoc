#ifndef ROBOTOC_RICCATI_RECURSION_HPP_
#define ROBOTOC_RICCATI_RECURSION_HPP_

#include "Eigen/Core"

#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/riccati/riccati_factorization.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/riccati_factorizer.hpp"
#include "robotoc/utils/aligned_vector.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/core/direction.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"
#include "robotoc/ocp/ocp_def.hpp"


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
  ~RiccatiRecursion() = default;

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
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] time_discretization Time discretization. 
  /// @param[in, out] kkt_matrix KKT matrix. 
  /// @param[in, out] kkt_residual KKT residual. 
  /// @param[in, out] factorization Riccati factorization. 
  ///
  void backwardRiccatiRecursion(const TimeDiscretization& time_discretization, 
                                KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual, 
                                RiccatiFactorization& factorization);

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] time_discretization Time discretization. 
  /// @param[in] kkt_matrix KKT matrix. 
  /// @param[in] kkt_residual KKT residual. 
  /// @param[in, out] factorization Riccati factorization. 
  /// @param[in] d Direction. 
  ///
  void forwardRiccatiRecursion(const TimeDiscretization& time_discretization, 
                               const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual, 
                               const RiccatiFactorization& factorization,
                               Direction& d) const;

  ///
  /// @brief Gets of the LQR policies over the horizon. 
  /// @return const reference to the LQR policies.
  ///
  const aligned_vector<LQRPolicy>& getLQRPolicy() const;

  ///
  /// @brief Resizes the internal data. 
  /// @param[in] time_discretization Time discretization. 
  ///
  void resizeData(const TimeDiscretization& time_discretization);

private:
  RiccatiFactorizer factorizer_;
  aligned_vector<LQRPolicy> lqr_policy_;
  aligned_vector<STOPolicy> sto_policy_;
  SplitRiccatiFactorization factorization_m_;

};

} // namespace robotoc

#endif // ROBOTkeps_OC_RICCATI_RECURSION_HPP_ 