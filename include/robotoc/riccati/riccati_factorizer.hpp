#ifndef ROBOTOC_RICCATI_FACTORIZER_HPP_ 
#define ROBOTOC_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/split_direction.hpp"
#include "robotoc/impulse/impulse_split_kkt_matrix.hpp"
#include "robotoc/impulse/impulse_split_kkt_residual.hpp"
#include "robotoc/impulse/impulse_split_direction.hpp"
#include "robotoc/ocp/switching_constraint_jacobian.hpp"
#include "robotoc/ocp/switching_constraint_residual.hpp"
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
  ///
  RiccatiFactorizer(const Robot& robot);

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
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization of the next stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in, out] riccati Riccati factorization of this impulse stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual, 
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs the forward Riccati recursion and computes the state 
  /// direction. 
  /// @param[in] kkt_matrix Split KKT matrix of this stage. 
  /// @param[in] kkt_residual Split KKT residual of this stage. 
  /// @param[in] lqr_policy LQR policy of this stage. 
  /// @param[in] d Split direction of this stage. 
  /// @param[in, out] d_next Split direction of the next stage. 
  /// @param[in] sto If true, the sensitivity w.r.t. the switching time is 
  /// considered. If false, it is not considered.
  ///
  template <typename SplitDirectionType>
  void forwardRiccatiRecursion(const SplitKKTMatrix& kkt_matrix, 
                               const SplitKKTResidual& kkt_residual,
                               const LQRPolicy& lqr_policy, SplitDirection& d, 
                               SplitDirectionType& d_next, 
                               const bool sto=false) const;

  ///
  /// @brief Performs the forward Riccati recursion and computes the state 
  /// direction. 
  /// @param[in] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in] d Split direction of this impulse stage. 
  /// @param[in, out] d_next Split direction of the next stage. 
  ///
  void forwardRiccatiRecursion(const ImpulseSplitKKTMatrix& kkt_matrix, 
                               const ImpulseSplitKKTResidual& kkt_residual,
                               const ImpulseSplitDirection& d, 
                               SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate. 
  /// @param[in] riccati Riccati factorization of this impulse stage. 
  /// @param[in, out] d Split direction. 
  /// @param[in] sto If true, the sensitivity w.r.t. the switching time is 
  /// considered. If false, it is not considered. Defalut is false.
  ///
  template <typename SplitDirectionType>
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      SplitDirectionType& d, 
                                      const bool sto=false);

  ///
  /// @brief Computes the Newton direction of the Lagrange multiplier with 
  /// respect to the switching constraint. 
  /// @param[in] c_riccati Riccati factorization for the switching constraint.
  /// @param[in, out] d Split direction of the this stage. 
  /// @param[in] sto If true, the sensitivity w.r.t. the switching time is 
  /// considered. If false, it is not considered. Defalut is false.
  ///
  static void computeLagrangeMultiplierDirection(
      const SplitConstrainedRiccatiFactorization& c_riccati,
      SplitDirection& d, const bool sto=false);

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::LLT<Eigen::MatrixXd> llt_, llt_s_;
  LQRPolicy lqr_policy_;
  BackwardRiccatiRecursionFactorizer backward_recursion_;

};

} // namespace robotoc

#include "robotoc/riccati/riccati_factorizer.hxx"

#endif // ROBOTOC_RICCATI_FACTORIZER_HPP_ 