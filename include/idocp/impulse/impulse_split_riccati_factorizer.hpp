#ifndef IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HPP_
#define IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/impulse/impulse_split_kkt_matrix.hpp"
#include "idocp/impulse/impulse_split_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hpp"


namespace idocp {

///
/// @class ImpulseSplitRiccatiFactorizer
/// @brief Riccati factorizer for ImpulseSplitOCP.
///
class ImpulseSplitRiccatiFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseSplitRiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitRiccatiFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitRiccatiFactorizer(
      const ImpulseSplitRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitRiccatiFactorizer& operator=(
      const ImpulseSplitRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitRiccatiFactorizer(
      ImpulseSplitRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitRiccatiFactorizer& operator=(
      ImpulseSplitRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization at the next time stage. 
  /// @param[in, out] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in, out] kkt_residual KKT residual at the current impulse stage. 
  /// @param[out] riccati Riccati factorization at the current impulse stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual, 
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs forward Riccati recursion and computes state direction. 
  /// @param[in] kkt_matrix KKT matrix at the current time stage. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in] d Split direction at the current time stage. 
  /// @param[out] d_next Split direction at the next time stage. 
  ///
  void forwardRiccatiRecursion(const ImpulseSplitKKTMatrix& kkt_matrix, 
                               const ImpulseSplitKKTResidual& kkt_residual,
                               const ImpulseSplitDirection& d, 
                               SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate vector. 
  /// @param[in] riccati Riccati factorization at the current impulse stage. 
  /// @param[out] d Split direction of the current impulse stage. 
  ///
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      ImpulseSplitDirection& d);

private:
  bool has_floating_base_;
  int dimv_;
  ImpulseBackwardRiccatiRecursionFactorizer backward_recursion_;

};

} // namespace idocp

#include "idocp/impulse/impulse_split_riccati_factorizer.hxx"

#endif // IDOCP_IMPULSE_SPLIT_RICCATI_FACTORIZER_HPP_ 