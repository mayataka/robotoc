#ifndef IDOCP_SPLIT_KKT_MATRIX_INVERTER_HPP_
#define IDOCP_SPLIT_KKT_MATRIX_INVERTER_HPP_

#include "Eigen/Dense"

#include "idocp/robot/robot.hpp"


namespace idocp {
///
/// @class SplitKKTMatrixInverter 
/// @brief Schur complement for SplitKKTMatrix.
///
class SplitKKTMatrixInverter {
public:
  ///
  /// @brief Construct a Schur complement.
  /// @param[in] robot Robot model. 
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

  ///
  /// @brief Computes the inverse of the split KKT matrix of the time stage. 
  /// @param[in] dt Time step of this time stage.
  /// @param[in] F Jacobian of the state equation.
  /// @param[in] Q Hessian of the Lagrangian.
  /// @param[in, out] KKT_mat_inv Inverse of the split KKT matrix.
  ///
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void invert(const double dt, const Eigen::MatrixBase<MatrixType1>& F, 
              const Eigen::MatrixBase<MatrixType2>& Q, 
              const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv);

  ///
  /// @brief Multiplies the Jacobian of the state equation to a matrix. 
  /// @param[in] dt Time step of this time stage.
  /// @param[in] F Jacobian of the state equation.
  /// @param[in] mat Multiplied matrix.
  /// @param[in, out] res Result.
  ///
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void multiplyF(const double dt, const Eigen::MatrixBase<MatrixType1>& F, 
                 const Eigen::MatrixBase<MatrixType2>& mat, 
                 const Eigen::MatrixBase<MatrixType3>& res);

  ///
  /// @brief Computes the inverse of the split KKT matrix of the time stage 
  /// including the switching constraint. 
  /// @param[in] dt Time step of this time stage.
  /// @param[in] F Jacobian of the state equation.
  /// @param[in] Pq Jacobian of the switching constraint.
  /// @param[in] Q Hessian of the Lagrangian.
  /// @param[in, out] KKT_mat_inv Inverse of the split KKT matrix.
  ///
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void invert(const double dt, const Eigen::MatrixBase<MatrixType1>& F, 
              const Eigen::MatrixBase<MatrixType2>& Pq, 
              const Eigen::MatrixBase<MatrixType3>& Q, 
              const Eigen::MatrixBase<MatrixType4>& KKT_mat_inv);

  ///
  /// @brief Multiplies the Jacobians of the state equation and the switching
  /// constraint to a matrix. 
  /// @param[in] dt Time step of this time stage.
  /// @param[in] F Jacobian of the state equation.
  /// @param[in] Pq Jacobian of the switching constraint.
  /// @param[in] mat Multiplied matrix.
  /// @param[in, out] res Result.
  ///
  template <typename MatrixType1, typename MatrixType2, typename MatrixType3, 
            typename MatrixType4>
  void multiplyFPq(const double dt, const Eigen::MatrixBase<MatrixType1>& F, 
                   const Eigen::MatrixBase<MatrixType2>& Pq, 
                   const Eigen::MatrixBase<MatrixType3>& mat, 
                   const Eigen::MatrixBase<MatrixType4>& res);

  ///
  /// @brief Sets the regularization. 
  /// @param[in] reg Regularization factor. Must be non-negative. Default is 0.
  ///
  void setRegularization(const double reg=0);

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