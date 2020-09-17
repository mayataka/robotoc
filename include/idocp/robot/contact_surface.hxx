#ifndef IDOCP_CONTACT_SURFACE_HXX_ 
#define IDOCP_CONTACT_SURFACE_HXX_ 

#include "idocp/robot/contact_surface.hpp"

namespace idocp {

inline ContactSurface::ContactSurface(const Eigen::Vector3d& origin, 
                                      const Eigen::Vector3d& normal_vector)
  : origin_(origin), 
    normal_vector_(normal_vector) {
}


inline ContactSurface::ContactSurface() {
  : origin_(origin), 
    normal_vector_(normal_vector) {
}


inline ContactSurface::~ContactSurface() {
}


inline bool ContactSurface::doesPenetrate(
    const Eigen::Vector3d& coordinate) const {
  if (point.coeff(2) < 0) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace idocp

#endif // IDOCP_CONTACT_SURFACE_HXX_ 