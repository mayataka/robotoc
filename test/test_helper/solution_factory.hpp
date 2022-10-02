#ifndef ROBOTOC_TEST_HELPER_SOLUTION_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_SOLUTION_FACTORY_HPP_

#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/time_discretization.hpp"
#include "robotoc/core/solution.hpp"


namespace robotoc {
namespace testhelper {

Solution CreateSolution(const Robot& robot, const int N);

Solution CreateSolution(const Robot& robot,  
                        const std::shared_ptr<ContactSequence>& contact_sequence,
                        const TimeDiscretization& time_discretization);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_SOLUTION_FACTORY_HPP_ 