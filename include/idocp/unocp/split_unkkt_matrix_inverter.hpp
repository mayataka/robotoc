#ifndef IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HPP_
#define IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitKKTMatrixInverter 
/// @brief Schur complement for SplitKKTMatrix.
///
class SplitUnKKTMatrixInverter {
public:
  ///
  /// @brief Construct a Schur complement.
  /// @param[in] robot Robot model. 
  ///
  SplitUnKKTMatrixInverter(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitUnKKTMatrixInverter();

  ///
  /// @brief Destructor. 
  ///
  ~SplitUnKKTMatrixInverter();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitUnKKTMatrixInverter(const SplitUnKKTMatrixInverter&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitUnKKTMatrixInverter& operator=(
      const SplitUnKKTMatrixInverter&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  SplitUnKKTMatrixInverter(SplitUnKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnKKTMatrixInverter& operator=(
      SplitUnKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Computes the inverse of the split KKT matrix of the time stage. 
  /// @param[in] dt Time step of this time stage.
  /// @param[in] Q Split KKT matrix.
  /// @param[in, out] KKT_mat_inv Inverse of the split KKT matrix.
  ///
  template <typename MatrixType1, typename MatrixType2>
  void invert(const double dt, const Eigen::MatrixBase<MatrixType1>& Q,
              const Eigen::MatrixBase<MatrixType2>& KKT_mat_inv);

private:
  Eigen::LLT<Eigen::MatrixXd> llt_Q_, llt_S_;
  Eigen::MatrixXd FQinv_, S_;
  int dimv_, dimx_, dimQ_, dimKKT_;

};

} // namespace idocp 

#include "idocp/unocp/split_unkkt_matrix_inverter.hxx"

#endif // IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HPP_ 