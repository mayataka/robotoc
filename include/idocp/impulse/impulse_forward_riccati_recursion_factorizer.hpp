#ifndef IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"


namespace idocp {

///
/// @class ImpulseForwardRiccatiRecursionFactorizer
/// @brief Factorizer of the backward Riccati recursion for SplitOCP.
///
class ImpulseForwardRiccatiRecursionFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseForwardRiccatiRecursionFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseForwardRiccatiRecursionFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseForwardRiccatiRecursionFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseForwardRiccatiRecursionFactorizer(
      const ImpulseForwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseForwardRiccatiRecursionFactorizer& operator=(
      const ImpulseForwardRiccatiRecursionFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseForwardRiccatiRecursionFactorizer(
      ImpulseForwardRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseForwardRiccatiRecursionFactorizer& operator=(
      ImpulseForwardRiccatiRecursionFactorizer&&) noexcept = default;

  ///
  /// @brief Factorize the KKT matrix and KKT residual for the backward Riccati  
  /// recursion, i.e., factorize matrices F, H, G and vectors lx, lu.
  /// @param[in] riccati_next Riccati factorization at the next time stage.
  /// @param[in] dtau Time step between the current time stage and the next 
  /// time stage.
  /// @param[in] kkt_matrix The KKT matrix.
  /// @param[in] kkt_residual The KKT residual.
  ///
  void factorizeStateTransition(const RiccatiFactorization& riccati, 
                                const ImpulseKKTMatrix& kkt_matrix, 
                                const ImpulseKKTResidual& kkt_residual, 
                                RiccatiFactorization& riccati_next);

  void factorizeStateConstraintFactorization(
      const RiccatiFactorization& riccati, const ImpulseKKTMatrix& kkt_matrix, 
      RiccatiFactorization& riccati_next);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd NApBKt_;

};

} // namespace idocp

#include "idocp/impulse/impulse_forward_riccati_recursion_factorizer.hxx"

#endif // IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 