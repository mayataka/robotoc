#ifndef ROBOTOC_TEST_HELPER_KKT_FACTORY_HPP_
#define ROBOTOC_TEST_HELPER_KKT_FACTORY_HPP_

#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/planner/contact_sequence.hpp"
#include "robotoc/ocp/ocp.hpp"
#include "robotoc/core/kkt_matrix.hpp"
#include "robotoc/core/kkt_residual.hpp"


namespace robotoc {
namespace testhelper {

SplitKKTMatrix CreateSplitKKTMatrix(const Robot& robot, const double dt);

SplitKKTMatrix CreateSplitKKTMatrix(const Robot& robot);

SplitKKTResidual CreateSplitKKTResidual(const Robot& robot);

KKTMatrix CreateKKTMatrix(const Robot& robot, const int N);

KKTMatrix CreateKKTMatrix(const Robot& robot, 
                          const std::shared_ptr<ContactSequence>& contact_sequence,
                          const TimeDiscretization& time_discretization);

KKTResidual CreateKKTResidual(const Robot& robot, const int N);

KKTResidual CreateKKTResidual(const Robot& robot, 
                              const std::shared_ptr<ContactSequence>& contact_sequence,
                              const TimeDiscretization& time_discretization);

} // namespace testhelper
} // namespace robotoc

#endif // ROBOTOC_TEST_HELPER_KKT_FACTORY_HPP_