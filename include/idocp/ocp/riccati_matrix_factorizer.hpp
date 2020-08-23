#ifndef IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_
#define IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_


#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class RiccatiMatrixFactorizer {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

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

  template <typename MatrixType1, typename MatrixType2>
  void setIntegrationSensitivities(
      const Eigen::MatrixBase<MatrixType1>& dintegrate_dq, 
      const Eigen::MatrixBase<MatrixType2>& dintegrate_dv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4, typename MatrixType5, typename MatrixType6, 
            typename MatrixType7, typename MatrixType8>
  void factorizeF(const double dtau, 
                  const Eigen::MatrixBase<MatrixType1>& Pqq_next,
                  const Eigen::MatrixBase<MatrixType2>& Pqv_next,
                  const Eigen::MatrixBase<MatrixType3>& Pvq_next,
                  const Eigen::MatrixBase<MatrixType4>& Pvv_next,
                  const Eigen::MatrixBase<MatrixType5>& Qqq,
                  const Eigen::MatrixBase<MatrixType6>& Qqv,
                  const Eigen::MatrixBase<MatrixType7>& Qvq,
                  const Eigen::MatrixBase<MatrixType8>& Qvv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void factorizeH(const double dtau, 
                  const Eigen::MatrixBase<MatrixType1>& Pqv_next,
                  const Eigen::MatrixBase<MatrixType2>& Pvv_next,
                  const Eigen::MatrixBase<MatrixType3>& Qqa,
                  const Eigen::MatrixBase<MatrixType4>& Qva);

  template <typename MatrixType1, typename MatrixType2>
  void factorizeG(const double dtau, 
                  const Eigen::MatrixBase<MatrixType1>& Pvv_next,
                  const Eigen::MatrixBase<MatrixType2>& Qaa);

private:
  bool has_floating_base_;
  int dimv_;
  Eigen::MatrixXd dintegrate_dq_, dintegrate_dv_;

};

} // namespace idocp

#include "idocp/ocp/riccati_matrix_factorizer.hxx"

#endif // IDOCP_RICCATI_MATRIX_FACTORIZER_HPP_