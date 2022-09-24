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


ImpulseStatus CreateActiveImpulseStatus(const Robot& robot, const double time_step) {
  auto impulse_status = robot.createImpulseStatus();
  for (int i=0; i<robot.contactFrames().size(); ++i) {
    impulse_status.setContactPlacement(i, Eigen::Vector3d::Random());
    impulse_status.setContactPlacement(i, SE3::Random());
  }
  impulse_status.setRandom();
  if (!impulse_status.hasActiveImpulse()) {
    impulse_status.activateImpulse(0);
  }
  return impulse_status;
};

} // namespace testhelper
} // namespace robotoc