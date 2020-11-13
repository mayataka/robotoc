#ifndef IDOCP_BACKWARD_RICCATI_RECURSION_HPP_ 
#define IDOCP_BACKWARD_RICCATI_RECURSION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/lqr_state_feedback_policy.hpp"


namespace idocp {

///
/// @class BackwardRiccatiRecursion
/// @brief Factorizer of the Riccati backward recursion for SplitOCP.
///
class BackwardRiccatiRecursion {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  BackwardRiccatiRecursion(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  BackwardRiccatiRecursion();

  ///
  /// @brief Destructor. 
  ///
  ~BackwardRiccatiRecursion();
 
  ///
  /// @brief Default copy constructor. 
  ///
  BackwardRiccatiRecursion(const BackwardRiccatiRecursion&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  BackwardRiccatiRecursion& operator=(const BackwardRiccatiRecursion&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  BackwardRiccatiRecursion(BackwardRiccatiRecursion&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  BackwardRiccatiRecursion& operator=(BackwardRiccatiRecursion&&) noexcept = default;

  void factorizeKKTMatrix(const RiccatiFactorization& riccati_next, 
                          const double dtau, KKTMatrix& kkt_matrix,  
                          KKTResidual& kkt_residual);

  void factorizeRiccatiFactorization(const RiccatiFactorization& riccati_next, 
                                     const KKTMatrix& kkt_matrix, 
                                     const KKTResidual& kkt_residual,
                                     const LQRStateFeedbackPolicy& lqr_policy,
                                     const double dtau, 
                                     RiccatiFactorization& riccati);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_, BtPq_, BtPv_, GK_;

};

} // namespace idocp

#include "idocp/ocp/backward_riccati_recursion.hxx"

#endif // IDOCP_BACKWARD_RICCATI_RECURSION_HPP_