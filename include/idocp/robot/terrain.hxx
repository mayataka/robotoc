#ifndef IDOCP_TERRAIN_HXX_
#define IDOCP_TERRAIN_HXX_

#include "idocp/robot/terrain.hpp"

namespace idocp {

inline Terrain::Terrain(const Eigen::Vector3d& origin, 
                        const Eigen::Vector3d& normal_vector)
  origin_(origin), 
  normal_vector_(normal_vector) {
}


inline Terrain::Terrain() {
}


inline Terrain::~Terrain() {
}


inline bool Terrain::doesPointPenetrate(const Eigen::Vector3d& point) const {
  if (point.coeff(2) < 0) {
    return true;
  }
  else {
    return false;
  }
}

} // namespace idocp

#endif // IDOCP_TERRAIN_HXX_