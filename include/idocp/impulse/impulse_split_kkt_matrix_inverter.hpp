#ifndef IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HPP_
#define IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HPP_

#include "Eigen/Dense"
#include "Eigen/LU"


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
  ImpulseSplitKKTMatrixInverter(const int dimv, const int max_dimf);

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
  void invert(const Eigen::MatrixBase<MatrixType1>& FC,
              const Eigen::MatrixBase<MatrixType2>& Q,
              const Eigen::MatrixBase<MatrixType3>& KKT_mat_inv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void multiplyFC(const Eigen::MatrixBase<MatrixType1>& FC, 
                  const Eigen::MatrixBase<MatrixType2>& mat, 
                  const Eigen::MatrixBase<MatrixType3>& res);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int dimv_, max_dimf_;
  bool has_floating_base_;
  Eigen::MatrixXd S_, FC_Qinv_;
  Eigen::LLT<Eigen::MatrixXd> llt_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_kkt_matrix_inverter.hxx"

#endif // IDOCP_IMPULSE_SPLIT_KKT_MATRIX_INVERTER_HPP_ 