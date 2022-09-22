#ifndef ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_

#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/core/direction.hpp"


namespace robotoc {
namespace testhelper {

Direction CreateDirection(const Robot& robot, const int N, 
                          const int max_num_impulse=0);

Direction CreateDirection(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence, 
                          const double T, const int N, 
                          const int max_num_impulse, const double t);


} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_ 