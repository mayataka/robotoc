#ifndef IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"


namespace idocp {

///
/// @class BackwardUnRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion for SplitOCP.
///
class BackwardUnRiccatiRecursionFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
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
  /// @brief Factorize the KKT matrix and KKT residual for the backward Riccati  
  /// recursion, i.e., factorize matrices F, H, G and vectors lx, lu.
  /// @param[in] riccati_next Riccati factorization at the next time stage.
  /// @param[in] dtau Time step between the current time stage and the next 
  /// time stage.
  /// @param[in, out] unkkt_matrix The condensed KKT matrix.
  /// @param[in, out] unkkt_residual The condensed KKT residual.
  ///
  void factorizeKKTMatrix(const SplitRiccatiFactorization& riccati_next, 
                          const double dtau, SplitUnKKTMatrix& unkkt_matrix,  
                          SplitUnKKTResidual& unkkt_residual);

  ///
  /// @brief Factorize the Riccati factorization matrix P and vector s.
  /// @param[in] riccati_next Riccati factorization at the next time stage.
  /// @param[in] unkkt_matrix The condensed KKT matrix factorized by 
  /// BackwardUnRiccatiRecursionFactorizer::factorizeRiccatiFactorization().
  /// @param[in] unkkt_residual The condensed KKT residual factorized by
  /// BackwardUnRiccatiRecursionFactorizer::factorizeRiccatiFactorization().
  /// @param[in] lqr_policy The state feedback control policy of the LQR 
  /// subproblem.
  /// @param[in] dtau Time step between the current time stage and the next 
  /// time stage.
  /// @param[out] riccati The Riccati factorization at this time stage.
  ///
  void factorizeRiccatiFactorization(
      const SplitRiccatiFactorization& riccati_next, 
      const SplitUnKKTMatrix& unkkt_matrix, 
      const SplitUnKKTResidual& unkkt_residual, 
      const LQRStateFeedbackPolicy& lqr_policy, const double dtau, 
      SplitRiccatiFactorization& riccati);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_;
  Eigen::MatrixXd GK_;

};

} // namespace idocp

#include "idocp/unocp/backward_unriccati_recursion_factorizer.hxx"

#endif // IDOCP_BACKWARD_UNRICCATI_RECURSION_FACTORIZER_HPP_ 