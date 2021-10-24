#ifndef ROBOTOC_LIE_DERIVATIVE_INVERTER_HPP_
#define ROBOTOC_LIE_DERIVATIVE_INVERTER_HPP_

#include "Eigen/Core"


namespace robotoc {

///
/// @class LieDerivativeInverter
/// @brief A class that computes the inverse of the Lie derivative of SE(3).
///
class LieDerivativeInverter {
public:
  ///
  /// @brief Default constructor.  
  ///
  LieDerivativeInverter();

  ///
  /// @brief Destructor. 
  ///
  ~LieDerivativeInverter();

  ///
  /// @brief Default copy constructor. 
  ///
  LieDerivativeInverter(const LieDerivativeInverter&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  LieDerivativeInverter& operator=(const LieDerivativeInverter&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  LieDerivativeInverter(LieDerivativeInverter&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  LieDerivativeInverter& operator=(LieDerivativeInverter&&) noexcept = default;

  ///
  /// @brief Computes the inverse of the Lie derivatives.
  /// @param[in] Lie_der The Lie derivative matrix. 
  /// Size must be larger than 6 x 6.
  /// @param[out] Lie_der_inv The inverse of the Lie derivatives. 
  /// Size must be larger than 6 x 6. Result is stored in the 6 x 6 top left 
  /// corner.
  ///
  template <typename MatrixType1, typename MatrixType2>
  void computeLieDerivativeInverse(
      const Eigen::MatrixBase<MatrixType1>& Lie_der,
      const Eigen::MatrixBase<MatrixType2>& Lie_der_inv);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::Matrix3d mat_3d_tmp_;

};

} // namespace robotoc

#include "robotoc/robot/lie_derivative_inverter.hxx"

#endif // ROBOTOC_LIE_DERIVATIVE_INVERTER_HPP_