#ifndef IDOCP_RAISIM_WRAPPER_HPP_
#define IDOCP_RAISIM_WRAPPER_HPP_

#include "Eigen/Core"

namespace idocp {

namespace raisimwrapper {

void pino2rai(const Eigen::VectorXd& q_pino, const Eigen::VectorXd& v_pino, 
              Eigen::VectorXd& q_raisim, Eigen::VectorXd& v_raisim);

void pino2rai(const Eigen::VectorXd& u_pinocchio, Eigen::VectorXd& u_raisim);

void rai2pino(const Eigen::VectorXd& q_raisim, const Eigen::VectorXd& v_raisim, 
              Eigen::VectorXd& q_pinocchio, Eigen::VectorXd& v_pinocchio);

} // namespace raisimwrapper
} // namespace idocp

#endif // IDOCP_RAISIM_WRAPPER_HPP_