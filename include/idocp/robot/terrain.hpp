#ifndef IDOCP_TERRAIN_HPP_
#define IDOCP_TERRAIN_HPP_

#include <vector>

#include "Eigen/Core"

namespace idocp {

class Terrain {
public:
  Terrain(const Eigen::Vector3d& origin=Eigen::Vector3d::Zero(), 
          const Eigen::Vector3d& normal_vector=Eigen::Vector3d::Zero());

  Terrain();
 
  ~Terrain();

  Terrain(const Terrain&) = default;

  Terrain& operator=(const Terrain&) = default;

  Terrain(Terrain&&) noexcept = default;

  Terrain& operator=(Terrain&&) noexcept = default;

  bool doesPointPenetrate(const Eigen::Vector3d& point) const;

private:
  Eigen::Vector3d origin_, normal_vector_;

};

} // namespace idocp

#include "idocp/robot/terrain.hxx"

#endif // IDOCP_TERRAIN_HPP_