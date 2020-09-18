#ifndef IDOCP_CONTACT_SURFACE_HPP_
#define IDOCP_CONTACT_SURFACE_HPP_

#include "Eigen/Core"

namespace idocp {

class ContactSurface {
public:
  ContactSurface(const Eigen::Vector3d& origin, 
                 const Eigen::Vector3d& normal_vector);

  ContactSurface();
 
  ~ContactSurface();

  ContactSurface(const ContactSurface&) = default;

  ContactSurface& operator=(const ContactSurface&) = default;

  ContactSurface(ContactSurface&&) noexcept = default;

  ContactSurface& operator=(ContactSurface&&) noexcept = default;

  bool doesPenetrate(const Eigen::Vector3d& coordinate) const;

  double normalContactForce(const Eigen::Vector3d& f) const;

  void setOrigin(const Eigen::Vector3d& origin);

  void setNormalVector(const Eigen::Vector3d& normal_vector);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::Vector3d origin_, normal_vector_;

};

} // namespace idocp

#include "idocp/robot/contact_surface.hxx"

#endif // IDOCP_CONTACT_SURFACE_HPP_ 