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
#include "idocp/ocp/backward_riccati_recursion.hpp"


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

  void setExistStateConstraint(const bool exist_state_constraint);

  void backwardRiccatiRecursion(const RiccatiFactorization& riccati_next, 
                                const double dtau, KKTMatrix& kkt_matrix, 
                                KKTResidual& kkt_residual,  
                                RiccatiFactorization& riccati);

  void forwardRiccatiRecursionParallel(KKTMatrix& kkt_matrix, 
                                       KKTResidual& kkt_residual);

  static void forwardRiccatiRecursionSerialInitial(
      const KKTMatrix& kkt_matrix, const KKTResidual& kkt_residual, 
      RiccatiFactorization& riccati);

  void forwardRiccatiRecursionSerial(const RiccatiFactorization& riccati, 
                                     const KKTMatrix& kkt_matrix, 
                                     const KKTResidual& kkt_residual, 
                                     RiccatiFactorization& riccati_next);

  template <typename SplitDirectionType>
  void forwardRiccatiRecursion(const KKTMatrix& kkt_matrix, 
                               const KKTResidual& kkt_residual,
                               const SplitDirection& d, const double dtau,
                               SplitDirectionType& d_next) const;

  template <typename VectorType>
  void computeStateDirection(const RiccatiFactorization& riccati, 
                             const Eigen::MatrixBase<VectorType>& dx0,
                             SplitDirection& d);

  void computeCostateDirection(const RiccatiFactorization& riccati, 
                               SplitDirection& d) const;

  void computeControlInputDirection(const RiccatiFactorization& riccati, 
                                    SplitDirection& d) const;

  template <typename MatrixType>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType>& K);

  template <typename MatrixType1, typename MatrixType2>
  void getStateFeedbackGain(const Eigen::MatrixBase<MatrixType1>& Kq,
                            const Eigen::MatrixBase<MatrixType2>& Kv);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_, exist_state_constraint_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::LLT<Eigen::MatrixXd> llt_;
  LQRStateFeedbackPolicy lqr_policy_;
  BackwardRiccatiRecursion backward_recursion_;
  Eigen::MatrixXd GinvBt_, BGinvBt_, NApBKt_;

};

} // namespace idocp

#include "idocp/ocp/riccati_factorizer.hxx"

#endif // IDOCP_RICCATI_FACTORIZER_HPP_