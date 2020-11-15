#ifndef IDOCP_IMPULSE_RICCATI_FACTORIZER_HPP_
#define IDOCP_IMPULSE_RICCATI_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"
#include "idocp/impulse/impulse_split_direction.hpp"
#include "idocp/impulse/impulse_backward_riccati_recursion_factorizer.hpp"
#include "idocp/impulse/impulse_forward_riccati_recursion_factorizer.hpp"


namespace idocp {

///
/// @class ImpulseRiccatiFactorizer
/// @brief Riccati factorizer for SplitImpulseOCP.
///
class ImpulseRiccatiFactorizer {
public:
  ///
  /// @brief Construct factorizer.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseRiccatiFactorizer(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseRiccatiFactorizer();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseRiccatiFactorizer();
 
  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseRiccatiFactorizer(const ImpulseRiccatiFactorizer&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseRiccatiFactorizer& operator=(const ImpulseRiccatiFactorizer&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseRiccatiFactorizer(ImpulseRiccatiFactorizer&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseRiccatiFactorizer& operator=(ImpulseRiccatiFactorizer&&) noexcept = default;

  void backwardRiccatiRecursion(const RiccatiFactorization& riccati_next, 
                                ImpulseKKTMatrix& kkt_matrix, 
                                ImpulseKKTResidual& kkt_residual, 
                                RiccatiFactorization& riccati);

  void forwardRiccatiRecursionSerial(const RiccatiFactorization& riccati, 
                                     const ImpulseKKTMatrix& kkt_matrix, 
                                     const ImpulseKKTResidual& kkt_residual, 
                                     RiccatiFactorization& riccati_next);

  template <typename MatrixType1, typename MatrixType2>
  void backwardStateConstraintFactorization(
      const ImpulseKKTMatrix& kkt_matrix, 
      const Eigen::MatrixBase<MatrixType1>& T_next,  
      const Eigen::MatrixBase<MatrixType2>& T) const;

  template <typename VectorType>
  static void computeStateDirection(
      const RiccatiFactorization& riccati, 
      const Eigen::MatrixBase<VectorType>& dx0, ImpulseSplitDirection& d);

  static void computeCostateDirection(const RiccatiFactorization& riccati, 
                                      ImpulseSplitDirection& d);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_;
  ImpulseBackwardRiccatiRecursionFactorizer backward_recursion_;
  ImpulseForwardRiccatiRecursionFactorizer forward_recursion_;

};

} // namespace idocp

#include "idocp/impulse/impulse_riccati_factorizer.hxx"

#endif // IDOCP_IMPULSE_RICCATI_FACTORIZER_HPP_ 