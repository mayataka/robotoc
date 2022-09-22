#include "robotoc/dynamics/switching_constraint_data.hpp"

#include <cassert>

namespace robotoc {

SwitchingConstraintData::SwitchingConstraintData(const Robot& robot)
  : q(Eigen::VectorXd::Zero(robot.dimq())),
    dq(Eigen::VectorXd::Zero(robot.dimv())),
    PqT_xi(Eigen::VectorXd::Zero(robot.dimv())) {
}


SwitchingConstraintData::SwitchingConstraintData()
  : q(),
    dq(),
    PqT_xi() {
}

} // namespace robotoc