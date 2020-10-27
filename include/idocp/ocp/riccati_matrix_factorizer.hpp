#ifndef IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_
#define IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_matrix.hpp"
#include "idocp/ocp/kkt_residual.hpp"
#include "idocp/ocp/riccati_factorization.hpp"
#include "idocp/ocp/riccati_gain.hpp"


namespace idocp {

class RiccatiMatrixFactorizer {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  RiccatiMatrixFactorizer(const Robot& robot);

  // Default constructor.
  RiccatiMatrixFactorizer();

  // Destructor.
  ~RiccatiMatrixFactorizer();
 
  // Use default copy constructor.
  RiccatiMatrixFactorizer(const RiccatiMatrixFactorizer&) = default;

  // Use default copy operator.
  RiccatiMatrixFactorizer& operator=(const RiccatiMatrixFactorizer&) = default;

  // Use default move constructor.
  RiccatiMatrixFactorizer(RiccatiMatrixFactorizer&&) noexcept = default;

  // Use default move assign operator.
  RiccatiMatrixFactorizer& operator=(RiccatiMatrixFactorizer&&) noexcept = default;

  void factorizeMatrices(const RiccatiFactorization& riccati_next, 
                         const double dtau, KKTMatrix& kkt_matrix, 
                         KKTResidual& kkt_residual);

  void factorizeRecursion(const RiccatiFactorization& riccati_next, 
                          const KKTMatrix& kkt_matrix, 
                          const KKTResidual& kkt_residual,
                          const RiccatiGain& gain, 
                          const double dtau,
                          RiccatiFactorization& riccati) const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  bool has_floating_base_;
  int dimv_, dimu_;
  static constexpr int kDimFloatingBase = 6;
  Eigen::MatrixXd AtPqq_, AtPqv_, AtPvq_, AtPvv_, BtPq_, BtPv_;
};

} // namespace idocp

#include "idocp/ocp/riccati_matrix_factorizer.hxx"

#endif // IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_