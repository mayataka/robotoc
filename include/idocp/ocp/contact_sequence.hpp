#ifndef IDOCP_CONTACT_SEQUENCE_HPP_
#define IDOCP_CONTACT_SEQUENCE_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {
class ContactSequence {
private:
  ContactSequence(const Robot& robot, const int N)
    : contact_sequence_(N, robot.is_each_contact_active()) {
  }

  ~ContactSequence() {
  }

  const std::vector<bool>& contact_status(const int i) const {
    return contact_sequence_[i];
  }

public:
  std::vector<std::vector<bool>> contact_sequence_;

};

} // namespace idocp 

#endif // IDOCP_CONTACT_SEQUENCE_HPP_