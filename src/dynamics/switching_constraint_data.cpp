#include "robotoc/dynamics/switching_constraint_data.hpp"

#include <cassert>

namespace robotoc {

SwitchingConstraintData::SwitchingConstraintData(const Robot& robot)
  : q(Eigen::VectorXd::Zero(robot.dimq())),
    dq(Eigen::VectorXd::Zero(robot.dimv())),
    PqT_xi(Eigen::VectorXd::Zero(robot.dimv())),
    Pq_full_(Eigen::MatrixXd::Zero(robot.max_dimf(), robot.dimv())),
    dimv_(robot.dimv()),
    dims_(0) {
}


SwitchingConstraintData::SwitchingConstraintData()
  : q(),
    dq(),
    PqT_xi(),
    Pq_full_(),
    dimv_(0),
    dims_(0) {
}

} // namespace robotoc