#ifndef IDOCP_TEST_HELPER_COST_FACTORY_HPP_
#define IDOCP_TEST_HELPER_COST_FACTORY_HPP_

#include <memory>

#include "idocp/cost/cost_function.hpp"


namespace idocp {
namespace testhelper {

std::shared_ptr<CostFunction> CreateCost(const Robot& robot);

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_COST_FACTORY_HPP_