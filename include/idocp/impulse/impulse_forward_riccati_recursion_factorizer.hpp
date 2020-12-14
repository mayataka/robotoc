#ifndef IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HPP_ 
#define IDOCP_IMPULSE_FORWARD_RICCATI_RECURSION_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"


namespace idocp {

///
/// @class ImpulseForwardRiccatiRecursionFactorizer
/// @brief Factorizer of the forward Riccati recursion for SplitOCP.
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
  /// @brief Factorize the state transition matrix and vector in terms of the 
  /// initial state.  
  /// @param[in] riccati Riccati factorization at the current impulse stage.
  /// @param[in] kkt_matrix The KKT matrix at the current impulse stage.
  /// @param[in] kkt_residual The KKT residual at the current impulse stage.
  /// @param[out] riccati_next Riccati factorization at the next time stage.
  ///
  void factorizeStateTransition(const RiccatiFactorization& riccati, 
                                const ImpulseSplitKKTMatrix& kkt_matrix, 
                                const ImpulseSplitKKTResidual& kkt_residual, 
                                RiccatiFactorization& riccati_next);

  ///
  /// @brief Factorize the pure-state constraint factorization matrix.  
  /// @param[in] riccati Riccati factorization at the current time stage.
  /// @param[in] kkt_matrix The KKT matrix the the current impulse stage.
  /// @param[out] riccati_next Riccati factorization at the next time stage.
  ///
  void factorizeStateConstraintFactorization(
      const RiccatiFactorization& riccati, 
      const ImpulseSplitKKTMatrix& kkt_matrix, 
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