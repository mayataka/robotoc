#ifndef IDOCP_SE3_HPP_
#define IDOCP_SE3_HPP_

#include <Eigen/Core>
#include "pinocchio/spatial/se3.hpp"


namespace idocp {

///
/// @typedef SE3
/// @brief Using pinocchio::SE3 without its namespace. 
///
using SE3 = pinocchio::SE3;

///
/// @brief Applies Log6 map that transforms the SE3 into 6-dimensional vector.
/// @param[in] SE3_obj SE3 object
/// @return Transformed 6-dimensional vector.
///
template <typename SE3Type>
const Eigen::Matrix<double, 6, 1> Log6Map(const SE3Type& SE3_obj) {
  return pinocchio::log6(SE3_obj).toVector();
}

///
/// @brief Computes the Jacobian of the Log6 map.
/// @param[in] SE3_obj SE3 object
/// @param[in, out] J Resultant Jacobian. Size must be 6x6.
///
template <typename SE3Type, typename MatrixType>
void computeJLog6Map(const SE3Type& SE3_obj, 
                     const Eigen::MatrixBase<MatrixType>& J) {
  assert(J.rows() == 6);
  assert(J.cols() == 6);
  pinocchio::Jlog6(SE3_obj, const_cast<Eigen::MatrixBase<MatrixType>&>(J));
}

} // namespace idocp

#endif // IDOCP_SE3_HPP_ 