#ifndef IDOCP_IMPULSE_RICCATI_MATRIX_FACTORIZER_HPP_
#define IDOCP_IMPULSE_RICCATI_MATRIX_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/impulse/impulse_kkt_matrix.hpp"
#include "idocp/impulse/impulse_kkt_residual.hpp"


namespace idocp {

class ImpulseRiccatiMatrixFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  ImpulseRiccatiMatrixFactorizer(const Robot& robot);

  // Default constructor.
  ImpulseRiccatiMatrixFactorizer();

  // Destructor.
  ~ImpulseRiccatiMatrixFactorizer();
 
  // Use default copy constructor.
  ImpulseRiccatiMatrixFactorizer(const ImpulseRiccatiMatrixFactorizer&) 
      = default;

  // Use default copy operator.
  ImpulseRiccatiMatrixFactorizer& operator=(const ImpulseRiccatiMatrixFactorizer&) 
      = default;

  // Use default move constructor.
  ImpulseRiccatiMatrixFactorizer(ImpulseRiccatiMatrixFactorizer&&) noexcept 
      = default;

  // Use default move assign operator.
  ImpulseRiccatiMatrixFactorizer& operator=(ImpulseRiccatiMatrixFactorizer&&) noexcept 
      = default;

  void factorize(const ImpulseKKTMatrix& kkt_matrix,
                 const ImpulseKKTResidual& kkt_residual,
                 const RiccatiFactorization& riccati_next,
                 RiccatiFactorization& riccati);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd ATPqq_, ATPqv_, ATPvq_, ATPvv_;

};

} // namespace idocp

#include "idocp/impulse/impulse_riccati_matrix_factorizer.hxx"

#endif // IDOCP_IMPULSE_RICCATI_MATRIX_FACTORIZER_HPP_ 