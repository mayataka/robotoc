#ifndef IDOCP_SPLIT_UNRICCATI_FACTORIZER_HPP_ 
#define IDOCP_SPLIT_UNRICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/unocp/split_unkkt_matrix.hpp"
#include "idocp/unocp/split_unkkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/unocp/backward_unriccati_recursion_factorizer.hpp"

#include <limits>
#include <cmath>


namespace idocp {

///
/// @class SplitUnRiccatiFactorizer
/// @brief Riccati factorizer for a time stage.
///
class SplitUnRiccatiFactorizer {
public:
  ///
  /// @brief Constructs a factorizer.
  /// @param[in] robot Robot model. 
  ///
  SplitUnRiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitUnRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnRiccatiFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  SplitUnRiccatiFactorizer(const SplitUnRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitUnRiccatiFactorizer& operator=(
      const SplitUnRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitUnRiccatiFactorizer(SplitUnRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnRiccatiFactorizer& operator=(
      SplitUnRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] dt Time step between the current time stage and the next 
  /// @param[in, out] unkkt_matrix KKT matrix at the this time stage. 
  /// @param[in, out] unkkt_residual KKT residual at the this time stage. 
  /// @param[out] riccati Riccati factorization at the this time stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                const double dt, 
                                SplitUnKKTMatrix& unkkt_matrix, 
                                SplitUnKKTResidual& unkkt_residual,  
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs forward Riccati recursion and computes state direction. 
  /// @param[in] unkkt_residual KKT residual at the current time stage. 
  /// @param[in] d Split direction at the current time stage. 
  /// @param[in] dt Time step between the current time stage and the next 
  /// @param[out] d_next Split direction at the next time stage. 
  ///
  void forwardRiccatiRecursion(const SplitUnKKTResidual& unkkt_residual,
                               SplitDirection& d, const double dt, 
                               SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate vector. 
  /// @param[in] riccati Riccati factorization at the current stage. 
  /// @param[in, out] d Split direction of the current this time stage. 
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      SplitDirection& d);

  ///
  /// @brief Getter of the state feedback gain of the LQR subproblem. 
  /// @param[in] K The state feedback gain. 
  ///
  template <typename MatrixType>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType>& K) const;

  ///
  /// @brief Getter of the state feedback gain of the LQR subproblem. 
  /// @param[in] Kq The state feedback gain with respect to the configuration. 
  /// @param[in] Kv The state feedback gain with respect to the velocity. 
  ///
  template <typename MatrixType1, typename MatrixType2>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Kq,
                            const Eigen::MatrixBase<MatrixType2>& Kv) const;

private:
  int dimv_;
  Eigen::LLT<Eigen::MatrixXd> llt_;
  LQRStateFeedbackPolicy lqr_policy_;
  BackwardUnRiccatiRecursionFactorizer backward_recursion_;

};

} // namespace idocp

#include "idocp/unocp/split_unriccati_factorizer.hxx"

#endif // IDOCP_SPLIT_UNRICCATI_FACTORIZER_HPP_ 