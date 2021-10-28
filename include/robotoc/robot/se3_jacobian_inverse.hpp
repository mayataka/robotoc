#ifndef ROBOTOC_SE3_JACOBIAN_INVERSE_HPP_
#define ROBOTOC_SE3_JACOBIAN_INVERSE_HPP_

#include "Eigen/Core"


namespace robotoc {

///
/// @class SE3JacobianInverse
/// @brief A class that computes the inverse of the Jacobian of SE(3).
///
class SE3JacobianInverse {
public:
  ///
  /// @brief Default constructor.  
  ///
  SE3JacobianInverse();

  ///
  /// @brief Destructor. 
  ///
  ~SE3JacobianInverse();

  ///
  /// @brief Default copy constructor. 
  ///
  SE3JacobianInverse(const SE3JacobianInverse&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SE3JacobianInverse& operator=(const SE3JacobianInverse&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SE3JacobianInverse(SE3JacobianInverse&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SE3JacobianInverse& operator=(SE3JacobianInverse&&) noexcept = default;

  ///
  /// @brief Computes the inverse of the Jacobian of SE3.
  /// @param[in] Jac The Lie derivative matrix. 
  /// Size must be larger than 6 x 6.
  /// @param[out] Jac_inv The inverse of the Lie derivatives. 
  /// Size must be larger than 6 x 6. Result is stored in the 6 x 6 top left 
  /// corner.
  ///
  template <typename MatrixType1, typename MatrixType2>
  void compute(const Eigen::MatrixBase<MatrixType1>& Jac,
               const Eigen::MatrixBase<MatrixType2>& Jac_inv);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::Matrix3d mat_3d_tmp_;

};

} // namespace robotoc

#include "robotoc/robot/se3_jacobian_inverse.hxx"

#endif // ROBOTOC_SE3_JACOBIAN_INVERSE_HPP_