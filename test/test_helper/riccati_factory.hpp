#ifndef IDOCP_TEST_RICCATI_FACTORY_HPP_
#define IDOCP_TEST_RICCATI_FACTORY_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/riccati/split_riccati_factorization.hpp"


namespace idocp {
namespace testhelper {

SplitRiccatiFactorization CreateSplitRiccatiFactorization(const Robot& robot);

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_RICCATI_FACTORY_HPP_ 