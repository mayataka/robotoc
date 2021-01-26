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
  void invert(const double dtau, const Eigen::MatrixBase<MatrixType1>& F, 
              const Eigen::MatrixBase<MatrixType2>& Q, 
              const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void multiplyF(const double dtau, const Eigen::MatrixBase<MatrixType1>& F, 
                 const Eigen::MatrixBase<MatrixType2>& mat, 
                 const Eigen::MatrixBase<MatrixType3>& res);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void invert(const double dtau, const Eigen::MatrixBase<MatrixType1>& F, 
              const Eigen::MatrixBase<MatrixType2>& Pq, 
              const Eigen::MatrixBase<MatrixType3>& Q, 
              const Eigen::MatrixBase<MatrixType4>& KKT_mat_inv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void multiplyFPq(const double dtau, const Eigen::MatrixBase<MatrixType1>& F, 
                   const Eigen::MatrixBase<MatrixType2>& Pq, 
                   const Eigen::MatrixBase<MatrixType3>& mat, 
                   const Eigen::MatrixBase<MatrixType4>& res);

  void enableRegularization(const double reg=1.0e-09);

  void disableRegularization();

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_, dimu_, dimx_, dimQ_, dimKKT_, dimf_;
  bool has_floating_base_, regularization_;
  double reg_;
  Eigen::LLT<Eigen::MatrixXd> llt_Q_, llt_F_, llt_FPq_;
  Eigen::MatrixXd S_full_, Jac_Qinv_full_;

  Eigen::Block<Eigen::MatrixXd> S();

  const Eigen::Block<const Eigen::MatrixXd> S() const;

  Eigen::Block<Eigen::MatrixXd> Jac_Qinv();

  const Eigen::Block<const Eigen::MatrixXd> Jac_Qinv() const;

};

} // namespace idocp 

#include "idocp/ocp/split_kkt_matrix_inverter.hxx"

#endif // IDOCP_SPLIT_KKT_MATRIX_INVERTER_HPP_ 