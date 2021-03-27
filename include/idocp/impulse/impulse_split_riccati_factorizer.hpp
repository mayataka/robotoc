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
/// @brief Riccati factorizer of an impulse stage.
///
class ImpulseSplitRiccatiFactorizer {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
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
  /// @param[in] riccati_next Riccati factorization of the next time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in, out] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in, out] riccati Riccati factorization of this impulse stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                ImpulseSplitKKTMatrix& kkt_matrix, 
                                ImpulseSplitKKTResidual& kkt_residual, 
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs the forward Riccati recursion and computes the state 
  /// direction. 
  /// @param[in] kkt_matrix Split KKT matrix of this impulse stage. 
  /// @param[in] kkt_residual Split KKT residual of this impulse stage. 
  /// @param[in] d Split direction of this impulse stage. 
  /// @param[in, out] d_next Split direction of the next time stage. 
  ///
  void forwardRiccatiRecursion(const ImpulseSplitKKTMatrix& kkt_matrix, 
                               const ImpulseSplitKKTResidual& kkt_residual,
                               const ImpulseSplitDirection& d, 
                               SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate. 
  /// @param[in] riccati Riccati factorization of this impulse stage. 
  /// @param[in, out] d Split direction of this impulse stage. 
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