#ifndef ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_

#include "robotoc/robot/robot.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/ocp/direction.hpp"


namespace robotoc {
namespace testhelper {

Direction CreateDirection(const Robot& robot, const int N, 
                          const int max_num_impulse=0);

Direction CreateDirection(const Robot& robot, const ContactSequence& contact_sequence, 
                          const double T, const int N, const int max_num_impulse, const double t);


} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_ 