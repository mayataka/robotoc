#ifndef ROBOTOC_LIE_DERIVATIVE_INVERTER_HXX_
#define ROBOTOC_LIE_DERIVATIVE_INVERTER_HXX_

#include "robotoc/robot/lie_derivative_inverter.hpp"

#include <cassert>


namespace robotoc {

inline LieDerivativeInverter::LieDerivativeInverter() 
  : mat_3d_tmp_(Eigen::Matrix3d::Zero()) {
}


inline LieDerivativeInverter::~LieDerivativeInverter() {
}


template <typename MatrixType1, typename MatrixType2>
inline void LieDerivativeInverter::computeLieDerivativeInverse(
    const Eigen::MatrixBase<MatrixType1>& Lie_der,
    const Eigen::MatrixBase<MatrixType2>& Lie_der_inv) {
  assert(Lie_der.rows() >= 6);
  assert(Lie_der.cols() >= 6);
  assert(Lie_der_inv.rows() >= 6);
  assert(Lie_der_inv.cols() >= 6);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Lie_der_inv).template topLeftCorner<3, 3>().noalias()
      = Lie_der.template topLeftCorner<3, 3>().inverse();
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Lie_der_inv).template block<3, 3>(3, 3).noalias()
      = Lie_der.template block<3, 3>(3, 3).inverse();
  mat_3d_tmp_.noalias() = Lie_der.template block<3, 3>(0, 3) 
                            * Lie_der_inv.template block<3, 3>(3, 3);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Lie_der_inv).template block<3, 3>(0, 3).noalias()
      = - Lie_der_inv.template topLeftCorner<3, 3>() * mat_3d_tmp_;
}

} // namespace robotoc

#endif // ROBOTOC_LIE_DERIVATIVE_INVERTER_HXX_