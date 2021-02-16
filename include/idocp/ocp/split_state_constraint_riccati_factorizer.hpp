#ifndef IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_ 
#define IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

///
/// @class SplitStateConstraintRiccatiFactorizer
/// @brief Riccati factorized matrix and vector for a state only equality 
/// constraint split into each impulse stages.
///
class SplitStateConstraintRiccatiFactorizer {
public:
  ///
  /// @brief Allocate Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of possible impulse. Must be 
  /// more than 1 and less than N. 
  ///
  SplitStateConstraintRiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitStateConstraintRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~SplitStateConstraintRiccatiFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  SplitStateConstraintRiccatiFactorizer(
      const SplitStateConstraintRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitStateConstraintRiccatiFactorizer& operator=(
      const SplitStateConstraintRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitStateConstraintRiccatiFactorizer(
      SplitStateConstraintRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitStateConstraintRiccatiFactorizer& operator=(
      SplitStateConstraintRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Set the dimension of the impulse. 
  /// @param[in] impulse_status Impulse status. 
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set the dimension of the impulse 0. 
  ///
  void resetImpulseStatus();

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

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  static constexpr double kMindtau = std::numeric_limits<double>::epsilon();
//   static constexpr double kMindtau
//       = std::sqrt(std::numeric_limits<double>::epsilon());
  Eigen::LLT<Eigen::MatrixXd> llt_;
  LQRStateFeedbackPolicy lqr_policy_;
  BackwardRiccatiRecursionFactorizer backward_recursion_;
  ForwardRiccatiRecursionFactorizer forward_recursion_;
  Eigen::MatrixXd GinvBt_, BGinvBt_;

};

} // namespace idocp

#include "idocp/ocp/split_state_constraint_riccati_factorizer.hxx"

#endif // IDOCP_SPLIT_STATE_CONSTRAINT_RICCATI_FACTORIZER_HPP_ 