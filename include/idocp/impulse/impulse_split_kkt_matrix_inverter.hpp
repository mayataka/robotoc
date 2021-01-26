#ifndef IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HPP_
#define IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HPP_

#include "Eigen/Dense"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"


namespace idocp {
///
/// @class ImpulseSplitKKTMatrixInverter 
/// @brief Schur complement for ImpulseSplitKKTMatrix.
///
class ImpulseSplitKKTMatrixInverter {
public:
  ///
  /// @brief Construct a Schur complement.
  ///
  ImpulseSplitKKTMatrixInverter(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitKKTMatrixInverter();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitKKTMatrixInverter();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitKKTMatrixInverter(const ImpulseSplitKKTMatrixInverter&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitKKTMatrixInverter& operator=(
      const ImpulseSplitKKTMatrixInverter&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitKKTMatrixInverter(
      ImpulseSplitKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitKKTMatrixInverter& operator=(
      ImpulseSplitKKTMatrixInverter&&) noexcept = default;

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void invert(const Eigen::MatrixBase<MatrixType1>& Jac,
              const Eigen::MatrixBase<MatrixType2>& Q,
              const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void multiplyJac(const Eigen::MatrixBase<MatrixType1>& Jac, 
                  const Eigen::MatrixBase<MatrixType2>& mat, 
                  const Eigen::MatrixBase<MatrixType3>& res);

  void enableRegularization(const double reg=1.0e-09);

  void disableRegularization();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_, max_dimf_, dimQ_;
  bool has_floating_base_, regularization_;
  double reg_;
  Eigen::MatrixXd S_full_, Jac_Qinv_full_;
  Eigen::LLT<Eigen::MatrixXd> llt_;

  Eigen::Block<Eigen::MatrixXd> S();

  const Eigen::Block<const Eigen::MatrixXd> S() const;

  Eigen::Block<Eigen::MatrixXd> Jac_Qinv();

  const Eigen::Block<const Eigen::MatrixXd> Jac_Qinv() const;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_kkt_matrix_inverter.hxx"

#endif // IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HPP_ 