#ifndef ROBOTOC_UNCONSTR_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_
#define ROBOTOC_UNCONSTR_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"
#include "robotoc/riccati/lqr_policy.hpp"


namespace robotoc {

///
/// @class UnconstrBackwardRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion of a time stage.
///
class UnconstrBackwardRiccatiRecursionFactorizer {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
  ///
  UnconstrBackwardRiccatiRecursionFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrBackwardRiccatiRecursionFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrBackwardRiccatiRecursionFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrBackwardRiccatiRecursionFactorizer(
      const UnconstrBackwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrBackwardRiccatiRecursionFactorizer& operator=(
      const UnconstrBackwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrBackwardRiccatiRecursionFactorizer(
      UnconstrBackwardRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrBackwardRiccatiRecursionFactorizer& operator=(
      UnconstrBackwardRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Factorizes the split KKT matrix and split KKT residual of a time 
  /// stage for the backward Riccati recursion.
  /// @param[in] riccati_next Riccati factorization of the next time stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void factorizeKKTMatrix(const SplitRiccatiFactorization& riccati_next, 
                          const double dt, SplitKKTMatrix& kkt_matrix,  
                          SplitKKTResidual& kkt_residual);

  ///
  /// @brief Factorizes the Riccati factorization matrix and vector.
  /// @param[in] riccati_next Riccati factorization of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in] lqr_policy The state feedback control policy of the LQR 
  /// subproblem.
  /// @param[in] dt Time step of this time stage.
  /// @param[out] riccati The Riccati factorization of this time stage.
  ///
  void factorizeRiccatiFactorization(
      const SplitRiccatiFactorization& riccati_next, 
      SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
      const LQRPolicy& lqr_policy, const double dt, 
      SplitRiccatiFactorization& riccati);

private:
  int dimv_;
  Eigen::MatrixXd GK_;

};

} // namespace robotoc

#include "robotoc/riccati/unconstr_backward_riccati_recursion_factorizer.hxx"

#endif // ROBOTOC_UNCONSTR_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_