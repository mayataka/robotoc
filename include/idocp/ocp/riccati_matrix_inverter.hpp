#ifndef IDOCP_RICCATI_MATRIX_INVERTER_HPP_
#define IDOCP_RICCATI_MATRIX_INVERTER_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class RiccatiMatrixInverter {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  RiccatiMatrixInverter(const Robot& robot);

  // Default constructor.
  RiccatiMatrixInverter();

  // Destructor.
  ~RiccatiMatrixInverter();
 
  // Use default copy constructor.
  RiccatiMatrixInverter(const RiccatiMatrixInverter&) = default;

  // Use default copy operator.
  RiccatiMatrixInverter& operator=(const RiccatiMatrixInverter&) = default;

  // Use default move constructor.
  RiccatiMatrixInverter(RiccatiMatrixInverter&&) noexcept = default;

  // Use default move assign operator.
  RiccatiMatrixInverter& operator=(RiccatiMatrixInverter&&) noexcept = default;

  void setContactStatus(const Robot& robot);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void invert(const Eigen::MatrixBase<MatrixType1>& G,
              const Eigen::MatrixBase<MatrixType2>& Caf,
              const Eigen::MatrixBase<MatrixType3>& G_inv);

  template <typename MatrixType1, typename MatrixType2>
  void invert(const Eigen::MatrixBase<MatrixType1>& G,
              const Eigen::MatrixBase<MatrixType2>& Caf);

  template <typename MatrixType>
  void getInverseMatrix(const Eigen::MatrixBase<MatrixType>& G_inv);

  template <typename MatrixType>
  void firstOrderCorrection(const double dtau, 
                            const Eigen::MatrixBase<MatrixType>& dPvv,
                            const Eigen::MatrixBase<MatrixType>& G_inv);
  
private:
  int dimv_, dimf_, dimc_, dimaf_;
  Eigen::MatrixXd G_inv_, Sc_, G_inv_Caf_trans_;
};

} // namespace idocp

#include "idocp/ocp/riccati_matrix_inverter.hxx"

#endif // IDOCP_RICCATI_MATRIX_INVERTER_HPP_ 