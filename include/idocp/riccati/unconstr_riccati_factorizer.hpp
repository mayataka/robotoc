#ifndef IDOCP_UNCONSTR_RICCATI_FACTORIZER_HPP_ 
#define IDOCP_UNCONSTR_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"
#include "idocp/riccati/lqr_policy.hpp"
#include "idocp/riccati/unconstr_backward_riccati_recursion_factorizer.hpp"

#include <limits>
#include <cmath>


namespace idocp {

///
/// @class UnconstrRiccatiFactorizer
/// @brief Riccati factorizer for a time stage.
///
class UnconstrRiccatiFactorizer {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
  ///
  UnconstrRiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrRiccatiFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrRiccatiFactorizer(const UnconstrRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrRiccatiFactorizer& operator=(
      const UnconstrRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  UnconstrRiccatiFactorizer(UnconstrRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrRiccatiFactorizer& operator=(
      UnconstrRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] dt Time step between the current time stage and the next 
  /// @param[in, out] kkt_matrix KKT matrix at the this time stage. 
  /// @param[in, out] kkt_residual KKT residual at the this time stage. 
  /// @param[out] riccati Riccati factorization at the this time stage. 
  /// @param[out] lqr_policy The state feedback control policy of the LQR 
  /// subproblem.
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                const double dt, SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual,  
                                SplitRiccatiFactorization& riccati,
                                LQRPolicy& lqr_policy);

  ///
  /// @brief Performs forward Riccati recursion and computes state direction. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in] dt Time step between the current time stage and the next 
  /// @param[in] lqr_policy The state feedback control policy of the LQR 
  /// subproblem.
  /// @param[in, out] d Split direction at the current time stage. 
  /// @param[out] d_next Split direction at the next time stage. 
  ///
  void forwardRiccatiRecursion(const SplitKKTResidual& kkt_residual,
                               const double dt, const LQRPolicy& lqr_policy,
                               SplitDirection& d,  SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate vector. 
  /// @param[in] riccati Riccati factorization at the current stage. 
  /// @param[in, out] d Split direction of the current this time stage. 
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      SplitDirection& d);

private:
  int dimv_;
  Eigen::LLT<Eigen::MatrixXd> llt_;
  UnconstrBackwardRiccatiRecursionFactorizer backward_recursion_;

};

} // namespace idocp

#include "idocp/riccati/unconstr_riccati_factorizer.hxx"

#endif // IDOCP_UNCONSTR_RICCATI_FACTORIZER_HPP_ 