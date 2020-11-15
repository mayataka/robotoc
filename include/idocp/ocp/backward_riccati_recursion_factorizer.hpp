#ifndef IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"


namespace idocp {

///
/// @class BackwardRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion for SplitOCP.
///
class BackwardRiccatiRecursionFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  BackwardRiccatiRecursionFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  BackwardRiccatiRecursionFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~BackwardRiccatiRecursionFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  BackwardRiccatiRecursionFactorizer(
      const BackwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  BackwardRiccatiRecursionFactorizer& operator=(
      const BackwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  BackwardRiccatiRecursionFactorizer(
      BackwardRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  BackwardRiccatiRecursionFactorizer& operator=(
      BackwardRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Factorize the KKT matrix and KKT residual for the backward Riccati  
  /// recursion, i.e., factorize matrices F, H, G and vectors lx, lu.
  /// @param[in] riccati_next Riccati factorization at the next time stage.
  /// @param[in] dtau Time step between the current time stage and the next 
  /// time stage.
  /// @param[in, out] kkt_matrix The KKT matrix.
  /// @param[in, out] kkt_residual The KKT residual.
  ///
  void factorizeKKTMatrix(const RiccatiFactorization& riccati_next, 
                          const double dtau, KKTMatrix& kkt_matrix,  
                          KKTResidual& kkt_residual);

  ///
  /// @brief Factorize the Riccati factorization matrix P and vector s.
  /// @param[in] riccati_next Riccati factorization at the next time stage.
  /// @param[in] kkt_matrix The KKT matrix factorized by 
  /// BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization().
  /// @param[in] kkt_residual The KKT residual factorized by
  /// BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization().
  /// @param[in] lqr_policy The state feedback control policy of the LQR 
  /// subproblem.
  /// @param[in] dtau Time step between the current time stage and the next 
  /// time stage.
  /// @param[out] riccati The Riccati factorization at this time stage.
  ///
  void factorizeRiccatiFactorization(const RiccatiFactorization& riccati_next, 
                                     const KKTMatrix& kkt_matrix, 
                                     const KKTResidual& kkt_residual,
                                     const LQRStateFeedbackPolicy& lqr_policy,
                                     const double dtau, 
                                     RiccatiFactorization& riccati);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_, BtPq_, BtPv_, GK_;

};

} // namespace idocp

#include "idocp/ocp/backward_riccati_recursion_factorizer.hxx"

#endif // IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 