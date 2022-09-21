#ifndef ROBOTOC_TEST_HELPER_SOLUTION_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_SOLUTION_FACTORY_HPP_

#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/hybrid/contact_sequence.hpp"
#include "robotoc/core/solution.hpp"


namespace robotoc {
namespace testhelper {

Solution CreateSolution(const Robot& robot, const int N, 
                        const int max_num_impulse=0);

Solution CreateSolution(const Robot& robot,  
                        const std::shared_ptr<ContactSequence>& contact_sequence, 
                        const double T, const int N, 
                        const int max_num_impulse, const double t);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_SOLUTION_FACTORY_HPP_ 