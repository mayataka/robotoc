#ifndef IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"


namespace idocp {

///
/// @class ImpulseBackwardRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion of an impulse stage.
///
class ImpulseBackwardRiccatiRecursionFactorizer {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
  ///
  ImpulseBackwardRiccatiRecursionFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseBackwardRiccatiRecursionFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseBackwardRiccatiRecursionFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseBackwardRiccatiRecursionFactorizer(
      const ImpulseBackwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseBackwardRiccatiRecursionFactorizer& operator=(
      const ImpulseBackwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseBackwardRiccatiRecursionFactorizer(
      ImpulseBackwardRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseBackwardRiccatiRecursionFactorizer& operator=(
      ImpulseBackwardRiccatiRecursionFactorizer&&) noexcept = default;

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
  int dimv_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_;

};

} // namespace idocp

#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hxx"

#endif // IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 