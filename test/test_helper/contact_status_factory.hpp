#ifndef ROBOTOC_TEST_HELPER_CONTACT_STATUS_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_CONTACT_STATUS_FACTORY_HPP_

#include "robotoc/robot/robot.hpp"


namespace robotoc {
namespace testhelper {

ContactStatus CreateActiveContactStatus(const Robot& robot, const double time_step);

ImpactStatus CreateActiveImpactStatus(const Robot& robot, const double time_step);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_CONTACT_STATUS_FACTORY_HPP_