#ifndef ROBOTOC_SE3_JACOBIAN_INVERSE_HXX_
#define ROBOTOC_SE3_JACOBIAN_INVERSE_HXX_

#include "robotoc/robot/se3_jacobian_inverse.hpp"

#include <cassert>


namespace robotoc {

inline SE3JacobianInverse::SE3JacobianInverse() 
  : mat_3d_tmp_(Eigen::Matrix3d::Zero()) {
}


template <typename MatrixType1, typename MatrixType2>
inline void SE3JacobianInverse::compute(
    const Eigen::MatrixBase<MatrixType1>& Jac, 
    const Eigen::MatrixBase<MatrixType2>& Jac_inv) {
  assert(Jac.rows() >= 6);
  assert(Jac.cols() >= 6);
  assert(Jac_inv.rows() >= 6);
  assert(Jac_inv.cols() >= 6);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Jac_inv).template topLeftCorner<3, 3>().noalias()
      = Jac.template topLeftCorner<3, 3>().inverse();
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Jac_inv).template block<3, 3>(3, 3).noalias()
      = Jac.template block<3, 3>(3, 3).inverse();
  mat_3d_tmp_.noalias() = Jac.template block<3, 3>(0, 3) 
                            * Jac_inv.template block<3, 3>(3, 3);
  const_cast<Eigen::MatrixBase<MatrixType2>&>(Jac_inv).template block<3, 3>(0, 3).noalias()
      = - Jac_inv.template topLeftCorner<3, 3>() * mat_3d_tmp_;
}

} // namespace robotoc

#endif // ROBOTOC_SE3_JACOBIAN_INVERSE_HXX_