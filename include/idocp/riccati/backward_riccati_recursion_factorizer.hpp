#ifndef IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"


namespace idocp {

///
/// @class BackwardRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion.
///
class BackwardRiccatiRecursionFactorizer {
public:
  using MatrixXdRowMajor 
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
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
  /// @brief Factorizes the split KKT matrix and split KKT residual of 
  /// this impulse stage for the backward Riccati recursion.
  /// @param[in] riccati_next Riccati factorization of the next time stage.
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage.
  ///
  void factorizeKKTMatrix(const SplitRiccatiFactorization& riccati_next, 
                          ImpulseSplitKKTMatrix& kkt_matrix);

  ///
  /// @brief Factorizes the Riccati factorization matrix and vector.
  /// @param[in] riccati_next Riccati factorization of the next time stage.
  /// @param[in] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in] kkt_residual Split KKT residual of this time stage.
  /// @param[in] lqr_policy The state feedback control policy of the LQR 
  /// subproblem.
  /// @param[in] dt Time step of this time stage.
  /// @param[out] riccati The Riccati factorization of this time stage.
  ///
  void factorizeRiccatiFactorization(
      const SplitRiccatiFactorization& riccati_next, 
      const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
      const LQRPolicy& lqr_policy, const double dt, 
      SplitRiccatiFactorization& riccati);

  ///
  /// @brief Factorizes the Riccati factorization matrix and vector.
  /// @param[in] riccati_next Riccati factorization of the next time stage.
  /// @param[in] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in] kkt_residual Split KKT residual of this impulse stage.
  /// ImpulseBackwardRiccatiRecursionFactorizer::factorizeKKTMatrix().
  /// @param[out] riccati The Riccati factorization of this impulse stage.
  ///
  void factorizeRiccatiFactorization(
      const SplitRiccatiFactorization& riccati_next, 
      const ImpulseSplitKKTMatrix& kkt_matrix, 
      const ImpulseSplitKKTResidual& kkt_residual, 
      SplitRiccatiFactorization& riccati);

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_, BtPq_, BtPv_, GK_, 
                  Pqq_tmp_, Pvv_tmp_;

};

} // namespace idocp

#include "idocp/riccati/backward_riccati_recursion_factorizer.hxx"

#endif // IDOCP_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 