#ifndef IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HPP_
#define IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HPP_

#include <vector>

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitUnKKTMatrixInverter
/// @brief The matrix inveter for SplitUnKKTMatrix.
///
class SplitUnKKTMatrixInverter {
public:
  ///
  /// @brief Construct a KKT matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
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
  SplitUnKKTMatrixInverter& operator=(const SplitUnKKTMatrixInverter&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  SplitUnKKTMatrixInverter(SplitUnKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitUnKKTMatrixInverter& operator=(SplitUnKKTMatrixInverter&&) noexcept = default;

  ///
  /// @brief Invert the split KKT matrix. 
  /// @param[in] dtau Discretization size of the optimal control problem.
  /// @param[in] Q Hessian of the Lagrangian with respect to the primal 
  /// variables.
  /// @param[out] KKT_matrix_inverse Inverse of the KKT matrix. Size must 
  /// be SplitUnKKTMatrix::dimKKT() x SplitUnKKTMatrix::dimKKT().
  ///
  template <typename MatrixType1, typename MatrixType2>
  void invert(const double dtau,  
              const Eigen::MatrixBase<MatrixType1>& Q,
              const Eigen::MatrixBase<MatrixType2>& KKT_matrix_inverse);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::LLT<Eigen::MatrixXd> llt_Q_, llt_S_;
  Eigen::MatrixXd FQinv_, S_;
  int dimv_, dimx_, dimQ_, dimKKT_;

};

} // namespace idocp 

#include "idocp/unocp/split_unkkt_matrix_inverter.hxx"

#endif // IDOCP_SPLIT_UNKKT_MATRIX_INVERTER_HPP_ 