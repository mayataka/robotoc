#ifndef IDOCP_SPLIT_RICCATI_FACTORIZER_HPP_ 
#define IDOCP_SPLIT_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/ocp/backward_riccati_recursion_factorizer.hpp"
#include "idocp/ocp/split_constrained_riccati_factorization.hpp"

#include <limits>
#include <cmath>


namespace idocp {

///
/// @class SplitRiccatiFactorizer
/// @brief Riccati factorizer for a time stage.
///
class SplitRiccatiFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  SplitRiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~SplitRiccatiFactorizer();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitRiccatiFactorizer(const SplitRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitRiccatiFactorizer& operator=(const SplitRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitRiccatiFactorizer(SplitRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitRiccatiFactorizer& operator=(SplitRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] dt Time step between the current time stage and the next 
  /// @param[in, out] kkt_matrix KKT matrix at the current stage. 
  /// @param[in, out] kkt_residual KKT residual at the current stage. 
  /// @param[out] riccati Riccati factorization at the current stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                const double dt, SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual,  
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs the backward Riccati recursion with the transformed 
  /// constraint. 
  /// @param[in] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] dt Time step between the current time stage and the next 
  /// @param[in, out] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in, out] kkt_residual KKT residual at the current impulse stage. 
  /// @param[in] jac The Jacobian of the transformed pure-state equality 
  /// constraint.
  /// @param[in, out] c_riccati Constrianed Riccati factorization at the current stage. 
  /// @param[out] riccati Riccati factorization at the current stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                const double dt, SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual,  
                                const SplitStateConstraintJacobian& jac,
                                SplitRiccatiFactorization& riccati,
                                SplitConstrainedRiccatiFactorization& c_riccati);

  ///
  /// @brief Performs forward Riccati recursion and computes state direction. 
  /// @param[in] kkt_matrix KKT matrix at the current time stage. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in] d Split direction at the current time stage. 
  /// @param[in] dt Time step between the current time stage and the next 
  /// @param[out] d_next Split direction at the next time stage. 
  ///
  template <typename SplitDirectionType>
  void forwardRiccatiRecursion(const SplitKKTMatrix& kkt_matrix, 
                               const SplitKKTResidual& kkt_residual,
                               const double dt, SplitDirection& d, 
                               SplitDirectionType& d_next) const;

  ///
  /// @brief Computes the Newton direction of the costate vector. 
  /// @param[in] riccati Riccati factorization at the current stage. 
  /// @param[in, out] d Split direction of the current impulse stage. 
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      SplitDirection& d);

  ///
  /// @brief Computes the Newton direction of the Lagrange multiplier with 
  /// respect to the switching constraint. 
  /// @param[in] c_riccati Constrained Riccati factorization at the current 
  /// stage. 
  /// @param[in, out] d Split direction of the current stage. 
  ///
  static void computeLagrangeMultiplierDirection(
      const SplitConstrainedRiccatiFactorization& c_riccati,
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
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::LLT<Eigen::MatrixXd> llt_, llt_s_;
  LQRStateFeedbackPolicy lqr_policy_;
  BackwardRiccatiRecursionFactorizer backward_recursion_;

};

} // namespace idocp

#include "idocp/ocp/split_riccati_factorizer.hxx"

#endif // IDOCP_SPLIT_RICCATI_FACTORIZER_HPP_ 