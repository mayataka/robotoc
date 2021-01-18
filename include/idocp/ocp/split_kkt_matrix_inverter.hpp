#ifndef IDOCP_SPLIT_KKT_MATRIX_INVERTER_HPP_
#define IDOCP_SPLIT_KKT_MATRIX_INVERTER_HPP_

#include "Eigen/Dense"

#include "idocp/robot/robot.hpp"


namespace idocp {
///
/// @class SplitKKTMatrixInverter
/// @brief The KKT matrix of a time stage.
///
class SplitKKTMatrixInverter {
public:
  ///
  /// @brief Construct a KKT matrix.
  ///
  SplitKKTMatrixInverter(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitKKTMatrixInverter();

  ///
  /// @brief Destructor. 
  ///
  ~SplitKKTMatrixInverter();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitKKTMatrixInverter(const SplitKKTMatrixInverter&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitKKTMatrixInverter& operator=(const SplitKKTMatrixInverter&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  SplitKKTMatrixInverter(SplitKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitKKTMatrixInverter& operator=(SplitKKTMatrixInverter&&) noexcept = default;

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void invert(const double dtau, const Eigen::MatrixBase<MatrixType1>& Jac, 
              const Eigen::MatrixBase<MatrixType2>& Q, 
              const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void multiplyJac(const double dtau, const Eigen::MatrixBase<MatrixType1>& Jac, 
                  const Eigen::MatrixBase<MatrixType2>& mat, 
                  const Eigen::MatrixBase<MatrixType3>& res);


  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_, dimu_, dimx_, dimQ_, dimKKT_;
  bool has_floating_base_;
  Eigen::LLT<Eigen::MatrixXd> llt_Q_, llt_F_;
  Eigen::MatrixXd S_, Jac_Qinv_;

};

} // namespace idocp 

#include "idocp/ocp/split_kkt_matrix_inverter.hxx"

#endif // IDOCP_SPLIT_KKT_MATRIX_INVERTER_HPP_ 