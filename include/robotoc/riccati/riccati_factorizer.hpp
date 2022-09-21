#ifndef ROBOTOC_RICCATI_FACTORIZER_HPP_ 
#define ROBOTOC_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "robotoc/robot/robot.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/switching_constraint_jacobian.hpp"
#include "robotoc/core/switching_constraint_residual.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"
#include "robotoc/riccati/sto_policy.hpp"
#include "robotoc/riccati/backward_riccati_recursion_factorizer.hpp"
#include "robotoc/riccati/split_constrained_riccati_factorization.hpp"

#include <limits>
#include <cmath>


namespace robotoc {

///
/// @class RiccatiFactorizer
/// @brief Riccati factorizer.
///
class RiccatiFactorizer {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
  /// @param[in] max_dts0 Maximum magnitude of the nominal direction of 
  /// the switching time. Used in a heuristic regularization on the dynamic 
  /// programming recursion. Must be positive. Default is 0.1.
  ///
  RiccatiFactorizer(const Robot& robot, const double max_dts0=0.1);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiFactorizer(const RiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  RiccatiFactorizer& operator=(const RiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiFactorizer(RiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiFactorizer& operator=(RiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Sets the regularization on the STO.
  /// @param[in] max_dts0 Maximum magnitude of the nominal direction of 
  /// the switching time. Used in a heuristic regularization on the dynamic 
  /// programming recursion. Must be positive. 
  ///
  void setRegularization(const double max_dts0);

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization of the next stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this stage. 
  /// @param[in, out] riccati Riccati factorization of this stage. 
  /// @param[in, out] lqr_policy LQR policy of this stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual,  
                                SplitRiccatiFactorization& riccati,
                                LQRPolicy& lqr_policy);

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization of the next stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this stage. 
  /// @param[in, out] riccati Riccati factorization of this stage. 
  /// @param[in, out] lqr_policy LQR policy of this stage. 
  /// @param[in] sto If true, the STO sensitivities are also considered. 
  /// @param[in] has_next_sto_phase Flag for wheather this phase has the next 
  /// phase involving the STO problem.
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual,  
                                SplitRiccatiFactorization& riccati,
                                LQRPolicy& lqr_policy, const bool sto,
                                const bool has_next_sto_phase);

  ///
  /// @brief Performs the backward Riccati recursion for the phase transition. 
  /// @param[in] riccati Riccati factorization of this stage. 
  /// @param[in, out] riccati_m Data for modified Riccati factorization. 
  /// @param[in, out] sto_policy STO policy. 
  /// @param[in] has_next_sto_phase Flag for wheather this phase has the next 
  /// phase involving the STO problem.
  ///
  void backwardRiccatiRecursionPhaseTransition(
      const SplitRiccatiFactorization& riccati, 
      SplitRiccatiFactorization& riccati_m, STOPolicy& sto_policy, 
      const bool has_next_sto_phase) const;

  ///
  /// @brief Performs the backward Riccati recursion with the switching 
  /// constraint. 
  /// @param[in] riccati_next Riccati factorization of the next stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this stage. 
  /// @param[in] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in] sc_residual Residual of the switching constraint. 
  /// @param[in, out] riccati Riccati factorization of this stage. 
  /// @param[in, out] c_riccati Riccati factorization for the switching 
  /// constraint. 
  /// @param[in, out] lqr_policy LQR policy of this stage. 
  ///
  void backwardRiccatiRecursion(
      const SplitRiccatiFactorization& riccati_next, 
      SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
      const SwitchingConstraintJacobian& sc_jacobian, 
      const SwitchingConstraintResidual& sc_residual, 
      SplitRiccatiFactorization& riccati, 
      SplitConstrainedRiccatiFactorization& c_riccati, LQRPolicy& lqr_policy);

  ///
  /// @brief Performs the backward Riccati recursion with the switching 
  /// constraint. 
  /// @param[in] riccati_next Riccati factorization of the next stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this stage. 
  /// @param[in] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in] sc_residual Residual of the switching constraint. 
  /// @param[in, out] riccati Riccati factorization of this stage. 
  /// @param[in, out] c_riccati Riccati factorization for the switching 
  /// constraint. 
  /// @param[in, out] lqr_policy LQR policy of this stage. 
  /// @param[in] sto If true, the STO sensitivities are also considered. 
  /// @param[in] has_next_sto_phase Flag for wheather this phase has the next 
  /// phase involving the STO problem.
  ///
  void backwardRiccatiRecursion(
      const SplitRiccatiFactorization& riccati_next, 
      SplitKKTMatrix& kkt_matrix, SplitKKTResidual& kkt_residual, 
      const SwitchingConstraintJacobian& sc_jacobian, 
      const SwitchingConstraintResidual& sc_residual, 
      SplitRiccatiFactorization& riccati, 
      SplitConstrainedRiccatiFactorization& c_riccati, LQRPolicy& lqr_policy,
      const bool sto, const bool has_next_sto_phase);

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization of the next stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in, out] riccati Riccati factorization of this impulse stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual, 
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization of the next stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in, out] riccati Riccati factorization of this impulse stage. 
  /// @param[in] sto If true, the STO sensitivities are also considered. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual, 
                                SplitRiccatiFactorization& riccati,
                                const bool sto);

  /// 
  /// @brief Computes the switching time direction. 
  /// @param[in, out] sto_policy STO policy. 
  /// @param[in] d Split direction. 
  /// @param[in] has_prev_sto_phase Flag for wheather this phase has the 
  /// previous phase involving the STO problem.
  /// 
  static void computeSwitchingTimeDirection(const STOPolicy& sto_policy,
                                            SplitDirection& d,
                                            const bool has_prev_sto_phase);

  ///
  /// @brief Performs the forward Riccati recursion and computes the state 
  /// direction. 
  /// @param[in] kkt_matrix Split KKT matrix of this stage. 
  /// @param[in] kkt_residual Split KKT residual of this stage. 
  /// @param[in] lqr_policy LQR policy of this stage. 
  /// @param[in] d Split direction of this stage. 
  /// @param[in, out] d_next Split direction of the next stage. 
  /// @param[in] sto If true, the STO sensitivities are also considered. 
  /// @param[in] has_next_sto_phase Flag for wheather this phase has the next 
  /// phase involving the STO problem.
  ///
  void forwardRiccatiRecursion(const SplitKKTMatrix& kkt_matrix, 
                               const SplitKKTResidual& kkt_residual,
                               const LQRPolicy& lqr_policy, SplitDirection& d, 
                               SplitDirection& d_next, const bool sto,
                               const bool has_next_sto_phase) const;

  ///
  /// @brief Performs the forward Riccati recursion and computes the state 
  /// direction. 
  /// @param[in] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in] d Split direction of this impulse stage. 
  /// @param[in, out] d_next Split direction of the next stage. 
  ///
  void forwardRiccatiRecursion(const SplitKKTMatrix& kkt_matrix, 
                               const SplitKKTResidual& kkt_residual,
                               const SplitDirection& d, 
                               SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate. 
  /// @param[in] riccati Riccati factorization of this impulse stage. 
  /// @param[in, out] d Split direction. 
  /// @param[in] sto If true, the STO sensitivities are also considered. 
  /// @param[in] has_next_sto_phase Flag for wheather this phase has the next 
  /// phase involving the STO problem.
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      SplitDirection& d, const bool sto, 
                                      const bool has_next_sto_phase);

  ///
  /// @brief Computes the Newton direction of the costate. 
  /// @param[in] riccati Riccati factorization of this impulse stage. 
  /// @param[in, out] d Impulse split direction. 
  /// @param[in] sto If true, the STO sensitivities are also considered. 
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      SplitDirection& d, const bool sto);

  ///
  /// @brief Computes the Newton direction of the Lagrange multiplier with 
  /// respect to the switching constraint. 
  /// @param[in] c_riccati Riccati factorization for the switching constraint.
  /// @param[in, out] d Split direction of the this stage. 
  /// @param[in] sto If true, the STO sensitivities are also considered. 
  /// @param[in] has_next_sto_phase Flag for wheather this phase has the next 
  /// phase involving the STO problem.
  ///
  static void computeLagrangeMultiplierDirection(
      const SplitConstrainedRiccatiFactorization& c_riccati, SplitDirection& d,
      const bool sto, const bool has_next_sto_phase);

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  double max_dts0_, eps_;
  Eigen::LLT<Eigen::MatrixXd> llt_, llt_s_;
  LQRPolicy lqr_policy_;
  BackwardRiccatiRecursionFactorizer backward_recursion_;

  template <typename SplitDirectionType>
  void forwardRiccatiRecursion_impl(
      const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
      const LQRPolicy& lqr_policy, SplitDirection& d, SplitDirectionType& d_next, 
      const bool sto, const bool has_next_sto_phase) const {
    d.du.noalias()  = lqr_policy.K * d.dx;
    d.du.noalias() += lqr_policy.k;
    if (sto) {
        d.du.noalias() += lqr_policy.T * (d.dts_next-d.dts);
        if (has_next_sto_phase) {
        d.du.noalias() -= lqr_policy.W * d.dts_next;
        }
    }
    d_next.dx = kkt_residual.Fx;
    d_next.dx.noalias()   += kkt_matrix.Fxx * d.dx;
    d_next.dv().noalias() += kkt_matrix.Fvu * d.du;
    if (sto) {
        d_next.dx.noalias() += kkt_matrix.fx * (d.dts_next-d.dts);
    }
    d_next.dts = d.dts;
    d_next.dts_next = d.dts_next;
  }

};

} // namespace robotoc

#endif // ROBOTOC_RICCATI_FACTORIZER_HPP_ 