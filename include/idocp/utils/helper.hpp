#ifndef IDOCP_HELPER_HPP_
#define IDOCP_HELPER_HPP_

#include <string>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {
namespace helper {

std::string GetCurrentWorkingDirectory();

Eigen::VectorXd pino2rai_q(const Eigen::VectorXd& q_pinocchio);

Eigen::VectorXd rai2pino_q(const Eigen::VectorXd& q_raisim);

Eigen::VectorXd pino2rai_v(Robot& robot, const Eigen::VectorXd& q_pinocchio, 
                           const Eigen::VectorXd& v_pinocchio);

Eigen::VectorXd rai2pino_v(Robot& robot, const Eigen::VectorXd& q_pinocchio, 
                           const Eigen::VectorXd& v_raisim);

Eigen::VectorXd pino2rai_u(const Eigen::VectorXd& u_pinocchio);

} // namespace helper
} // namespace idocp 

#endif // IDOCP_HELPER_HPP_ 