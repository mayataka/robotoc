#include "contact_status_factory.hpp"

#include <string>

#include "robotoc/robot/robot.hpp"


namespace robotoc {
namespace testhelper {

ContactStatus CreateActiveContactStatus(const Robot& robot, const double time_step) {
  auto contact_status = robot.createContactStatus();
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    contact_status.setContactPlacement(i, Eigen::Vector3d::Random());
    contact_status.setContactPlacement(i, SE3::Random());
  }
  contact_status.setRandom();
  if (!contact_status.hasActiveContacts()) {
    contact_status.activateContact(0);
  }
  return contact_status;
};


ImpactStatus CreateActiveImpactStatus(const Robot& robot, const double time_step) {
  auto impact_status = robot.createImpactStatus();
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    impact_status.setContactPlacement(i, Eigen::Vector3d::Random());
    impact_status.setContactPlacement(i, SE3::Random());
  }
  impact_status.setRandom();
  if (!impact_status.hasActiveImpact()) {
    impact_status.activateImpact(0);
  }
  return impact_status;
};

} // namespace testhelper
} // namespace robotoc