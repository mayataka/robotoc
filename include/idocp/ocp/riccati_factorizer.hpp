#ifndef IDOCP_RICCATI_FACTORIZER_HPP_
#define IDOCP_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"
#include "idocp/ocp/backward_riccati_recursion_factorizer.hpp"
#include "idocp/ocp/forward_riccati_recursion_factorizer.hpp"


namespace idocp {

///
/// @class RiccatiFactorizer
/// @brief Riccati factorizer for SplitOCP.
///
class RiccatiFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  RiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  RiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~RiccatiFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  RiccatiFactorizer(const RiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  RiccatiFactorizer& operator=(const RiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  RiccatiFactorizer(RiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  RiccatiFactorizer& operator=(RiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Performs the backward Riccati recursion. 
  /// @param[in] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] dtau Time step between the current time stage and the next 
  /// @param[in, out] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in, out] kkt_residual KKT residual at the current impulse stage. 
  /// @param[out] riccati Riccati factorization at the current impulse stage. 
  ///
  void backwardRiccatiRecursion(const RiccatiFactorization& riccati_next, 
                                const double dtau, KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual,  
                                RiccatiFactorization& riccati);

  ///
  /// @brief Performs the parallel part of the forward Riccati recursion with 
  /// pure-state equality constraints. 
  /// @param[in, out] kkt_matrix KKT matrix at the current impulse stage. 
  /// @param[in, out] kkt_residual KKT residual at the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. Default is false.
  ///
  void forwardRiccatiRecursionParallel(KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual,
                                       const bool exist_state_constraint=false);

  ///
  /// @brief Serial part of the forward Riccati recursion with state 
  /// constraint. 
  ///
  static void forwardRiccatiRecursionSerialInitial(
      const RiccatiFactorization& riccati);

  ///
  /// @brief Performs the serial part of the forward Riccati recursion with 
  /// pure-state equality constraints. 
  /// @param[in] riccati Riccati factorization at the current time stage. 
  /// @param[in] kkt_matrix KKT matrix at the current time stage. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in] dtau Time step between the current time stage and the next 
  /// @param[out] riccati_next Riccati factorization at the next time stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. Default is false.
  ///
  void forwardRiccatiRecursionSerial(const RiccatiFactorization& riccati, 
                                     const KKTMatrix& kkt_matrix, 
                                     const KKTResidual& kkt_residual, 
                                     const double dtau,
                                     RiccatiFactorization& riccati_next,
                                     const bool exist_state_constraint=false);

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
      const KKTMatrix& kkt_matrix, const Eigen::MatrixBase<MatrixType1>& T_next,  
      const double dtau, const Eigen::MatrixBase<MatrixType2>& T) const;

  ///
  /// @brief This is unconstrained version of the forward Riccati recursion. 
  /// @param[in] kkt_matrix KKT matrix at the current time stage. 
  /// @param[in] kkt_residual KKT residual at the current time stage. 
  /// @param[in, out] d Split direction at the current time stage. 
  /// @param[in] dtau Time step between the current time stage and the next 
  /// @param[out] d_next Split direction at the next time stage. 
  ///
  void forwardRiccatiRecursion(const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual,
                               SplitDirection& d, const double dtau,
                               SplitDirection& d_next) const;

  ///
  /// @brief Computes the Newton direction of the state vector. 
  /// @param[in] riccati Riccati factorization at the current impulse stage. 
  /// @param[in] dx0 Direction of the state at the initial time stage. 
  /// @param[out] d Split direction of the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. Default is false.
  ///
  template <typename VectorType>
  static void computeStateDirection(const RiccatiFactorization& riccati, 
                                    const Eigen::MatrixBase<VectorType>& dx0,
                                    SplitDirection& d,
                                    const bool exist_state_constraint=false);

  ///
  /// @brief Computes the Newton direction of the costate vector. 
  /// @param[in] riccati Riccati factorization at the current impulse stage. 
  /// @param[in, out] d Split direction of the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. Default is false.
  ///
  static void computeCostateDirection(const RiccatiFactorization& riccati, 
                                      SplitDirection& d,
                                      const bool exist_state_constraint=false);

  ///
  /// @brief Computes the Newton direction of the control input vector. 
  /// @param[in] riccati Riccati factorization at the current impulse stage. 
  /// @param[in, out] d Split direction of the current impulse stage. 
  /// @param[in] exist_state_constraint If true, the factorization for state
  /// constraints are also performed. Default is false.
  ///
  void computeControlInputDirection(
      const RiccatiFactorization& riccati, SplitDirection& d,
      const bool exist_state_constraint=false) const;

  ///
  /// @brief Getter of the state feedback gain of the LQR subproblem. 
  /// @param[in] K The state feedback gain. 
  ///
  template <typename MatrixType>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType>& K);

  ///
  /// @brief Getter of the state feedback gain of the LQR subproblem. 
  /// @param[in] Kq The state feedback gain with respect to the configuration. 
  /// @param[in] Kv The state feedback gain with respect to the velocity. 
  ///
  template <typename MatrixType1, typename MatrixType2>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Kq,
                            const Eigen::MatrixBase<MatrixType2>& Kv);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_, exist_state_constraint;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::LLT<Eigen::MatrixXd> llt_;
  LQRStateFeedbackPolicy lqr_policy_;
  BackwardRiccatiRecursionFactorizer backward_recursion_;
  ForwardRiccatiRecursionFactorizer forward_recursion_;
  Eigen::MatrixXd GinvBt_, BGinvBt_;

};

} // namespace idocp

#include "idocp/ocp/riccati_factorizer.hxx"

#endif // IDOCP_RICCATI_FACTORIZER_HPP_