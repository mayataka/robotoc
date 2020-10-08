#ifndef IDOCP_CONTACT_SURFACE_HXX_ 
#define IDOCP_CONTACT_SURFACE_HXX_ 

#include "idocp/robot/contact_surface.hpp"

#include <stdexcept>
#include <iostream>

namespace idocp {

inline ContactSurface::ContactSurface(const Eigen::Vector3d& origin, 
                                      const Eigen::Vector3d& normal_vector)
  : origin_(origin), 
    normal_vector_(normal_vector) {
  try {
    if (normal_vector.squaredNorm() != 1.0) {
      throw std::out_of_range(
          "invalid argment: norm of normal vector must be 1");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
}


inline ContactSurface::ContactSurface()
  : origin_(Eigen::Vector3d::Zero()), 
    normal_vector_(Eigen::Vector3d::Zero()) {
  normal_vector_.coeffRef(2) = 1;
}


inline ContactSurface::~ContactSurface() {
}


inline bool ContactSurface::doesPenetrate(
    const Eigen::Vector3d& coordinate) const {
  if (coordinate.coeff(2) < 0) {
    return true;
  }
  else {
    return false;
  }
}


inline double ContactSurface::normalContactForce(
    const Eigen::Vector3d& f) const {
  return f.dot(normal_vector_);
}


inline void ContactSurface::setOrigin(const Eigen::Vector3d& origin) {
  origin_ = origin;
}


inline void ContactSurface::setNormalVector(
    const Eigen::Vector3d& normal_vector) {
  try {
    if (normal_vector.squaredNorm() != 1.0) {
      throw std::out_of_range(
          "invalid argment: norm of normal vector must be 1");
    }
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << '\n';
    std::exit(EXIT_FAILURE);
  }
  normal_vector_ = normal_vector;
}

} // namespace idocp

#endif // IDOCP_CONTACT_SURFACE_HXX_ 