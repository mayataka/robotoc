#ifndef IDOCP_TEST_HELPER_KKT_FACTORY_HPP_
#define IDOCP_TEST_HELPER_KKT_FACTORY_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"
#include "idocp/hybrid/hybrid_container.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"


namespace idocp {
namespace testhelper {

SplitKKTMatrix CreateSplitKKTMatrix(const Robot& robot, const double dt);

SplitKKTResidual CreateSplitKKTResidual(const Robot& robot);

SplitKKTResidual CreateSplitKKTResidual(const Robot& robot, 
                                        const ImpulseStatus& impulse_status);

KKTMatrix CreateKKTMatrix(const Robot& robot, const ContactSequence& contact_sequence, 
                          const int N, const int max_num_impulse);

KKTResidual CreateKKTResidual(const Robot& robot, const ContactSequence& contact_sequence, 
                              const int N, const int max_num_impulse);

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_KKT_FACTORY_HPP_