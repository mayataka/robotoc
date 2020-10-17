#ifndef IDOCP_SCHUR_COMPLEMENT_HPP_
#define IDOCP_SCHUR_COMPLEMENT_HPP_

#include "Eigen/Dense"


namespace idocp {
///
/// @class SchurComplement
/// @brief The KKT matrix of a time stage.
///
class SchurComplement {
public:
  ///
  /// @brief Construct a KKT matrix.
  ///
  SchurComplement(const int max_dimA, const int max_dimD);

  ///
  /// @brief Default constructor. 
  ///
  SchurComplement();

  ///
  /// @brief Destructor. 
  ///
  ~SchurComplement();

  ///
  /// @brief Default copy constructor. 
  ///
  SchurComplement(const SchurComplement&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SchurComplement& operator=(const SchurComplement&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  SchurComplement(SchurComplement&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SchurComplement& operator=(SchurComplement&&) noexcept = default;

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void invertWithZeroBottomRightCorner(const int dimA, const int dimD, 
                                       const Eigen::MatrixBase<MatrixType1>& A, 
                                       const Eigen::MatrixBase<MatrixType2>& C, 
                                       const Eigen::MatrixBase<MatrixType3>& Minv);

  template <typename MatrixType1, typename MatrixType2, typename MatrixType3>
  void invertWithZeroTopLeftCorner(const int dimA, const int dimD,
                                   const Eigen::MatrixBase<MatrixType1>& B,
                                   const Eigen::MatrixBase<MatrixType2>& D,
                                   const Eigen::MatrixBase<MatrixType3>& Minv);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  int max_dimA_, max_dimD_;
  Eigen::MatrixXd SA_, SD_, CAinv_, BDinv_;

};

} // namespace idocp 

#include "idocp/ocp/schur_complement.hxx"

#endif // IDOCP_SCHUR_COMPLEMENT_HPP_ 