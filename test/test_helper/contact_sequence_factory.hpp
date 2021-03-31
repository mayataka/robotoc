#ifndef IDOCP_TEST_HELPER_CONTACT_SEQUENCE_FACTORY_HPP_
#define IDOCP_TEST_HELPER_CONTACT_SEQUENCE_FACTORY_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/hybrid/contact_sequence.hpp"


namespace idocp {
namespace testhelper {

ContactSequence CreateContactSequence(const Robot& robot, const int N, 
                                      const int max_num_impulse,
                                      const double t0,
                                      const double event_period);

} // namespace testhelper
} // namespace idocp

#endif // IDOCP_TEST_HELPER_CONTACT_SEQUENCE_FACTORY_HPP_ 