#ifndef ROBOTOC_ROTATION_HPP_
#define ROBOTOC_ROTATION_HPP_

#include "Eigen/Core"
#include "Eigen/Geometry"

namespace robotoc {
namespace rotation {

///
/// @brief Convert quaternion vector (x, y, z, w) to a Rotation matrix.
/// @param[in] quat_xyzw Quaternion vector (x, y, z, w).
/// @return Rotation matrix.
///
template <typename VectorType>
inline Eigen::Matrix3d RotationMatrixFromQuaternion(
    const Eigen::MatrixBase<VectorType>& quat_xyzw) {
  assert(quat_xyzw.size() == 4);
  return Eigen::Quaterniond(quat_xyzw).toRotationMatrix();
}


///
/// @brief Convert a normal vector to its surface Rotation matrix.
/// @param[in] normal_vector Normal vector.
/// @return Rotation matrixo of the surface.
///
template <typename VectorType>
inline Eigen::Matrix3d RotationMatrixFromNormalVector(
    const Eigen::MatrixBase<VectorType>& normal_vector) {
  assert(normal_vector.size() == 3);
  Eigen::Matrix3d R;
  const double nx = normal_vector.coeff(0);
  const double ny = normal_vector.coeff(1);
  const double nz = normal_vector.coeff(2);
  const double nxny_norm = std::sqrt(nx*nx + ny*ny);
  constexpr double eps = std::numeric_limits<double>::epsilon();
  if (nxny_norm < eps) return Eigen::Matrix3d::Identity();

  R <<    ny/nxny_norm,   -nx/nxny_norm,         0.,
       nx*nz/nxny_norm, ny*nz/nxny_norm, -nxny_norm,
                    nx,              ny,         nz;
  return R;
}

///
/// @brief Convert a rotation matrix to a quaternion vector (x, y, z, w).
/// @param[in] R rotation matrix.
/// @return Quaternion vector (x, y, z, w).
///
template <typename MatrixType>
inline Eigen::Vector4d QuaternionFromRotationMatrix(
    const Eigen::MatrixBase<MatrixType>& R) {
  assert(R.rows() == 3);
  assert(R.cols() == 3);
  return Eigen::Quaterniond(R).coeffs();
}

///
/// @brief Convert a normal vector to its surface quaternion (x, y, z, w).
/// @param[in] normal_vector Normal vector.
/// @return Quaternion vector (x, y, z, w).
///
template <typename VectorType>
inline Eigen::Vector4d QuaternionFromNormalVector(
    const Eigen::MatrixBase<VectorType>& normal_vector) {
  return QuaternionFromRotationMatrix(RotationMatrixFromNormalVector(normal_vector));
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
inline void ProjectRotationMatrix(Eigen::Matrix3d& R, const ProjectionAxis axis) {
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