#ifndef ROBOTOC_TEST_HELPER_CONSTRAINTS_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_CONSTRAINTS_FACTORY_HPP_

#include <memory>

#include "robotoc/constraints/constraints.hpp"


namespace robotoc {
namespace testhelper {

std::shared_ptr<Constraints> CreateConstraints(const Robot& robot);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_CONSTRAINTS_FACTORY_HPP_