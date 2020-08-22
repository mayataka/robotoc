#ifndef IDOCP_CONTACT_SEQUENCE_HPP_
#define IDOCP_CONTACT_SEQUENCE_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {
class ContactSequence {
private:
  ContactSequence(const Robot& robot, const int N) {

  }

  ~ContactSequence();
};

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_HPP_