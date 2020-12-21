#ifndef IDOCP_SPLIT_RICCATI_FACTORIZER_HPP_ 
#define IDOCP_SPLIT_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/ocp/backward_riccati_recursion_factorizer.hpp"
#include "idocp/ocp/forward_riccati_recursion_factorizer.hpp"

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
  /// @param[in] dtau Time step between the current time stage and the next 
  /// @param[in, out] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in, out] kkt_residual KKT residual at the current impulse stage. 
  /// @param[out] riccati Riccati factorization at the current impulse stage. 
  ///
  void backwardRiccatiRecursion(const SplitRiccatiFactorization& riccati_next, 
                                const double dtau, SplitKKTMatrix& kkt_matrix, 
                                SplitKKTResidual& kkt_residual,  
                                SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs the parallel part of the forward Riccati recursion with 
  /// pure-state equality constraints. 
  /// @param[in, out] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in, out] kkt_residual KKT residual at the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. 
  ///
  void forwardRiccatiRecursionParallel(SplitKKTMatrix& kkt_matrix, 
                                       SplitKKTResidual& kkt_residual,
                                       const bool exist_state_constraint);

  ///
  /// @brief Checks the initial factorization of the serial part of the forward 
  /// Riccati recursion with pure-state equality constraints. 
  /// @param[in] riccati Riccati factorization at the current time stage. 
  ///
  static void forwardStateConstraintFactorizationInitial(
      const SplitRiccatiFactorization& riccati);

  ///
  /// @brief Performs the serial part of the forward Riccati recursion due to
  /// pure-state equality constraints. 
  /// @param[in] riccati Riccati factorization at the current time stage. 
  /// @param[in] kkt_matrix KKT matrix at the current time stage. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in] dtau Time step between the current time stage and the next 
  /// @param[out] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. 
  ///
  void forwardStateConstraintFactorization(
      const SplitRiccatiFactorization& riccati, 
      const SplitKKTMatrix& kkt_matrix, const SplitKKTResidual& kkt_residual, 
      const double dtau, SplitRiccatiFactorization& riccati_next, 
      const bool exist_state_constraint);

  ///
  /// @brief Performs the backward factorization of matrices related to the 
  /// pure-state equality constraints. 
  /// @param[in] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in] T_next A factorization at the next time stage. 
  /// @param[in] dtau Time step between the current time stage and the next 
  /// @param[out] T A factorization at the current impulse stage. 
  ///
  template <typename MatrixType1, typename MatrixType2>
  void backwardStateConstraintFactorization(
      const Eigen::MatrixBase<MatrixType1>& T_next, 
      const SplitKKTMatrix& kkt_matrix, const double dtau, 
      const Eigen::MatrixBase<MatrixType2>& T) const;

  ///
  /// @brief Performs forward Riccati recursion and computes state direction. 
  /// @param[in] kkt_matrix KKT matrix at the current time stage. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in] riccati_next Riccati factorization at the next stage. 
  /// @param[in] d Split direction at the current time stage. 
  /// @param[in] dtau Time step between the current time stage and the next 
  /// @param[out] d_next Split direction at the next time stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. 
  ///
  template <typename SplitDirectionType>
  void forwardRiccatiRecursion(const SplitKKTMatrix& kkt_matrix, 
                               const SplitKKTResidual& kkt_residual,
                               const SplitRiccatiFactorization& riccati_next,
                               const SplitDirection& d, 
                               const double dtau, 
                               SplitDirectionType& d_next, 
                               const bool exist_state_constraint) const;

  ///
  /// @brief Computes the Newton direction of the costate vector. 
  /// @param[in] riccati Riccati factorization at the current stage. 
  /// @param[in, out] d Split direction of the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed.
  ///
  static void computeCostateDirection(const SplitRiccatiFactorization& riccati, 
                                      SplitDirection& d,
                                      const bool exist_state_constraint);

  ///
  /// @brief Computes the Newton direction of the control input vector. 
  /// @param[in] riccati_next Riccati factorization at the next stage. 
  /// @param[in, out] d Split direction of the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. 
  ///
  void computeControlInputDirection(
      const SplitRiccatiFactorization& riccati_next, SplitDirection& d,
      const bool exist_state_constraint) const;

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

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_, is_dtau_sufficiently_positive_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  static constexpr double kMindtau
      = std::sqrt(std::numeric_limits<double>::epsilon());
  Eigen::LLT<Eigen::MatrixXd> llt_;
  LQRStateFeedbackPolicy lqr_policy_;
  BackwardRiccatiRecursionFactorizer backward_recursion_;
  ForwardRiccatiRecursionFactorizer forward_recursion_;
  Eigen::MatrixXd GinvBt_, BGinvBt_;

};

} // namespace idocp

#include "idocp/ocp/split_riccati_factorizer.hxx"

#endif // IDOCP_SPLIT_RICCATI_FACTORIZER_HPP_ 