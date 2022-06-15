#ifndef ROBOTOC_ROTATION_HPP_
#define ROBOTOC_ROTATION_HPP_

#include "Eigen/Core"
#include "Eigen/Geometry"

namespace robotoc {
namespace rotation {

///
/// @brief Convert quaternion vector (x, y, z, w) to Eigen::Quaterniond.
/// @param[in] quat_xyzw Quaternion vector (x, y, z, w).
/// @return Quaternion expressed as Eigen::Quaterniond.
///
template <typename VectorType>
inline Eigen::Quaterniond toQuaternion(const Eigen::MatrixBase<VectorType>& quat_xyzw) {
  assert(quat_xyzw.size() == 4);
  return Eigen::Quaterniond(quat_xyzw.coeff(3),  // w
                            quat_xyzw.coeff(0),  // x
                            quat_xyzw.coeff(1),  // y
                            quat_xyzw.coeff(2)); // z
}


///
/// @brief Convert quaternion vector (x, y, z, w) to a Rotation matrix.
/// @param[in] quat_xyzw Quaternion vector (x, y, z, w).
/// @return Rotation matrix.
///
template <typename VectorType>
inline Eigen::Matrix3d toRotationMatrix(const Eigen::MatrixBase<VectorType>& quat_xyzw) {
  assert(quat_xyzw.size() == 4);
  return toQuaternion(quat_xyzw).toRotationMatrix();
}


///
/// @enum ProjectionAxis
/// @brief Projection axis of the rotation. 
///
enum class ProjectionAxis {
  X, Y, Z,
};


///
/// @brief Projects a rotation matrix onto a specified axis.
/// @param[in, out] R Rotation matrix.
/// @param[in] axis Axis of the projection.
///
inline void projectRotationMatrix(Eigen::Matrix3d& R, const ProjectionAxis axis) {
  if (axis == ProjectionAxis::X) {
    const double norm = R.coeff(1, 1) * R.coeff(1, 1) + R.coeff(1, 2) + R.coeff(1, 2);
    R.array() /= norm;
    R.coeffRef(0, 0) = 1.0;
    R.coeffRef(0, 1) = 0.0;
    R.coeffRef(0, 2) = 0.0;
    R.coeffRef(1, 0) = 0.0;
    R.coeffRef(2, 0) = 0.0;
  }
  else if (axis == ProjectionAxis::Y) {
    const double norm = R.coeff(0, 0) * R.coeff(0, 0) + R.coeff(0, 2) + R.coeff(0, 2);
    R.array() /= norm;
    R.coeffRef(0, 1) = 0.0;
    R.coeffRef(1, 0) = 0.0;
    R.coeffRef(1, 1) = 1.0;
    R.coeffRef(1, 2) = 0.0;
    R.coeffRef(2, 1) = 0.0;
  }
  else { // ProjectionAxis::Z
    const double norm = R.coeff(0, 0) * R.coeff(0, 0) + R.coeff(0, 1) + R.coeff(0, 1);
    R.array() /= norm;
    R.coeffRef(0, 2) = 0.0;
    R.coeffRef(1, 2) = 0.0;
    R.coeffRef(2, 0) = 0.0;
    R.coeffRef(2, 1) = 0.0;
    R.coeffRef(2, 2) = 1.0;
  }
}

} // namespace rotation
} // namespace robotoc 

#endif // ROBOTOC_ROTATION_HPP_