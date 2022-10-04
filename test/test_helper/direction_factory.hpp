#ifndef ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_

#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/core/direction.hpp"


namespace robotoc {
namespace testhelper {

Direction CreateDirection(const Robot& robot, const int N);

Direction CreateDirection(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence,
                          const TimeDiscretization& time_discretization);


} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_DIRECTION_FACTORY_HPP_ 