#ifndef IDOCP_TEST_HELPER_CONSTRAINTS_FACTORY_HPP_
#define IDOCP_TEST_HELPER_CONSTRAINTS_FACTORY_HPP_

#include <memory>

#include "idocp/constraints/constraints.hpp"


namespace idocp {
namespace testhelper {

std::shared_ptr<Constraints> CreateConstraints(const Robot& robot);

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_CONSTRAINTS_FACTORY_HPP_