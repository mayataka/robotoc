#ifndef IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"


namespace idocp {

///
/// @class BackwardUnRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion of a time stage.
///
class BackwardUnRiccatiRecursionFactorizer {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
  ///
  BackwardUnRiccatiRecursionFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  BackwardUnRiccatiRecursionFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~BackwardUnRiccatiRecursionFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  BackwardUnRiccatiRecursionFactorizer(
      const BackwardUnRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  BackwardUnRiccatiRecursionFactorizer& operator=(
      const BackwardUnRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  BackwardUnRiccatiRecursionFactorizer(
      BackwardUnRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  BackwardUnRiccatiRecursionFactorizer& operator=(
      BackwardUnRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Factorizes the split KKT matrix and split KKT residual of a time 
  /// stage for the backward Riccati recursion.
  /// @param[in] riccati_next Riccati factorization of the next time stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in, out] unkkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] unkkt_residual Split KKT residual of this time stage.
  ///
  void factorizeKKTMatrix(const SplitRiccatiFactorization& riccati_next, 
                          const double dt, SplitUnKKTMatrix& unkkt_matrix,  
                          SplitUnKKTResidual& unkkt_residual);

  ///
  /// @brief Factorizes the Riccati factorization matrix and vector.
  /// @param[in] riccati_next Riccati factorization of the next time stage.
  /// @param[in] unkkt_matrix Split KKT matrix of this time stage.
  /// @param[in] unkkt_residual Split KKT residual of this time stage.
  /// @param[in] lqr_policy The state feedback control policy of the LQR 
  /// subproblem.
  /// @param[in] dt Time step of this time stage.
  /// @param[out] riccati The Riccati factorization of this time stage.
  ///
  void factorizeRiccatiFactorization(
      const SplitRiccatiFactorization& riccati_next, 
      const SplitUnKKTMatrix& unkkt_matrix, 
      const SplitUnKKTResidual& unkkt_residual, 
      const LQRPolicy& lqr_policy, const double dt, 
      SplitRiccatiFactorization& riccati);

private:
  int dimv_;
  Eigen::MatrixXd GK_;

};

} // namespace idocp

#include "idocp/unocp/backward_unriccati_recursion_factorizer.hxx"

#endif // IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HPP_ 