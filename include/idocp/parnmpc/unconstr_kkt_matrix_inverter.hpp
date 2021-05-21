#ifndef IDOCP_UNCONSTR_KKT_MATRIX_INVERTER_HPP_
#define IDOCP_UNCONSTR_KKT_MATRIX_INVERTER_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class UnconstrKKTMatrixInverter 
/// @brief Schur complement for SplitKKTMatrix for UnconstrParNMPC.
///
class UnconstrKKTMatrixInverter {
public:
  ///
  /// @brief Construct a Schur complement.
  /// @param[in] robot Robot model. 
  ///
  UnconstrKKTMatrixInverter(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  UnconstrKKTMatrixInverter();

  ///
  /// @brief Destructor. 
  ///
  ~UnconstrKKTMatrixInverter();

  ///
  /// @brief Default copy constructor. 
  ///
  UnconstrKKTMatrixInverter(const UnconstrKKTMatrixInverter&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  UnconstrKKTMatrixInverter& operator=(
      const UnconstrKKTMatrixInverter&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  UnconstrKKTMatrixInverter(UnconstrKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  UnconstrKKTMatrixInverter& operator=(
      UnconstrKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Computes the inverse of the split KKT matrix of the time stage. 
  /// @param[in] dt Time step of this time stage.
  /// @param[in] H Hessian of the KKT matrix.
  /// @param[in, out] KKT_mat_inv Inverse of the split KKT matrix.
  ///
  template <typename MatrixType1, typename MatrixType2>
  void invert(const double dt, const Eigen::MatrixBase<MatrixType1>& H,
              const Eigen::MatrixBase<MatrixType2>& KKT_mat_inv);

private:
  Eigen::LLT<Eigen::MatrixXd> llt_H_, llt_S_;
  Eigen::MatrixXd FHinv_, S_;
  int dimv_, dimx_, dimH_, dimkkt_;

};

} // namespace idocp 

#include "idocp/parnmpc/unconstr_kkt_matrix_inverter.hxx"

#endif // IDOCP_UNCONSTR_KKT_MATRIX_INVERTER_HPP_ 