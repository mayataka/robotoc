#ifndef IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"


namespace idocp {

///
/// @class ImpulseBackwardRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion for SplitImpulseOCP.
///
class ImpulseBackwardRiccatiRecursionFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
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
  /// @brief Factorize the KKT matrix and KKT residual for the backward Riccati  
  /// recursion, i.e., factorize matrices F, and vectors lx.
  /// @param[in] riccati_next Riccati factorization at the next time stage.
  /// @param[in, out] kkt_matrix The KKT matrix.
  ///
  void factorizeKKTMatrix(const RiccatiFactorization& riccati_next, 
                          ImpulseSplitKKTMatrix& kkt_matrix);

  ///
  /// @brief Factorize the Riccati factorization matrix P and vector s.
  /// @param[in] riccati_next Riccati factorization at the next time stage.
  /// @param[in] kkt_matrix The KKT matrix factorized by 
  /// BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization().
  /// @param[in] kkt_residual The KKT residual factorized by
  /// BackwardRiccatiRecursionFactorizer::factorizeRiccatiFactorization().
  /// @param[out] riccati The Riccati factorization at this time stage.
  ///
  void factorizeRiccatiFactorization(
      const RiccatiFactorization& riccati_next, 
      const ImpulseSplitKKTMatrix& kkt_matrix, 
      const ImpulseSplitKKTResidual& kkt_residual, 
      RiccatiFactorization& riccati);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_;

};

} // namespace idocp

#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hxx"

#endif // IDOCP_IMPULSE_BACKWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 