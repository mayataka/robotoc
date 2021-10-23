#ifndef ROBOTOC_TEST_RICCATI_FACTORY_HPP_
#define ROBOTOC_TEST_RICCATI_FACTORY_HPP_

#include "robotoc/robot/robot.hpp"
#include "robotoc/riccati/split_riccati_factorization.hpp"


namespace robotoc {
namespace testhelper {

SplitRiccatiFactorization CreateSplitRiccatiFactorization(const Robot& robot);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_RICCATI_FACTORY_HPP_ 