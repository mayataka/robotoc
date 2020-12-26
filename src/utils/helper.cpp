#include "idocp/utils/helper.hpp"

#include <iostream>
#include <unistd.h>


namespace idocp {
namespace helper {

std::string GetCurrentWorkingDirectory() {
  constexpr int _PATH_MAX = 255;
  char cwd[_PATH_MAX];
  if (getcwd(cwd, _PATH_MAX) == NULL) {
    std::cerr << "Failed in getting CWD!!" << std::endl;
  }
  return std::string(cwd, strlen(cwd));
}


Eigen::VectorXd pino2rai_q(const Eigen::VectorXd& q_pinocchio) {
  if (q_pinocchio.size() == 19) { // quadruped configuration
    Eigen::VectorXd q_raisim(Eigen::VectorXd::Zero(19));
    // Base frame position
    q_raisim.coeffRef(0)  = q_pinocchio.coeff(0);
    q_raisim.coeffRef(1)  = q_pinocchio.coeff(1);
    q_raisim.coeffRef(2)  = q_pinocchio.coeff(2);
    // Change quaternion order
    q_raisim.coeffRef(3)  = q_pinocchio.coeff(6);
    q_raisim.coeffRef(4)  = q_pinocchio.coeff(3);
    q_raisim.coeffRef(5)  = q_pinocchio.coeff(4);
    q_raisim.coeffRef(6)  = q_pinocchio.coeff(5);
    // LH 
    q_raisim.coeffRef(7)  = q_pinocchio.coeff(7);
    q_raisim.coeffRef(8)  = q_pinocchio.coeff(8);
    q_raisim.coeffRef(9)  = q_pinocchio.coeff(9);
    // Change leg order LH <-> RF
    q_raisim.coeffRef(10) = q_pinocchio.coeff(13);
    q_raisim.coeffRef(11) = q_pinocchio.coeff(14);
    q_raisim.coeffRef(12) = q_pinocchio.coeff(15);
    // Change leg order LH <-> RF
    q_raisim.coeffRef(13) = q_pinocchio.coeff(10);
    q_raisim.coeffRef(14) = q_pinocchio.coeff(11);
    q_raisim.coeffRef(15) = q_pinocchio.coeff(12);
    // RH
    q_raisim.coeffRef(16) = q_pinocchio.coeff(16);
    q_raisim.coeffRef(17) = q_pinocchio.coeff(17);
    q_raisim.coeffRef(18) = q_pinocchio.coeff(18);
    return q_raisim;
  }
  else {
    return q_pinocchio;
  }
}


Eigen::VectorXd rai2pino_q(const Eigen::VectorXd& q_raisim) {
  if (q_raisim.size() == 19) { // quadruped configuration
    Eigen::VectorXd q_pinocchio(Eigen::VectorXd::Zero(19));
    // Base frame position
    q_pinocchio.coeffRef(0)  = q_raisim.coeff(0);
    q_pinocchio.coeffRef(1)  = q_raisim.coeff(1);
    q_pinocchio.coeffRef(2)  = q_raisim.coeff(2);
    // Change quaternion order
    q_pinocchio.coeffRef(6)  = q_raisim.coeff(3);
    q_pinocchio.coeffRef(3)  = q_raisim.coeff(4);
    q_pinocchio.coeffRef(4)  = q_raisim.coeff(5);
    q_pinocchio.coeffRef(5)  = q_raisim.coeff(6);
    // LH 
    q_pinocchio.coeffRef(7)  = q_raisim.coeff(7);
    q_pinocchio.coeffRef(8)  = q_raisim.coeff(8);
    q_pinocchio.coeffRef(9)  = q_raisim.coeff(9);
    // Change leg order LH <-> RF
    q_pinocchio.coeffRef(13) = q_raisim.coeff(10);
    q_pinocchio.coeffRef(14) = q_raisim.coeff(11);
    q_pinocchio.coeffRef(15) = q_raisim.coeff(12);
    // Change leg order LH <-> RF
    q_pinocchio.coeffRef(10) = q_raisim.coeff(13);
    q_pinocchio.coeffRef(11) = q_raisim.coeff(14);
    q_pinocchio.coeffRef(12) = q_raisim.coeff(15);
    // RH
    q_pinocchio.coeffRef(16) = q_raisim.coeff(16);
    q_pinocchio.coeffRef(17) = q_raisim.coeff(17);
    q_pinocchio.coeffRef(18) = q_raisim.coeff(18);
    return q_pinocchio;
  }
  else {
    return q_raisim;
  }
}


Eigen::VectorXd pino2rai_v(Robot& robot, const Eigen::VectorXd& q_pinocchio, 
                           const Eigen::VectorXd& v_pinocchio) {
  if (v_pinocchio.size() == 18) { // quadruped velocity
    assert(q_pinocchio.size() === 19);
    robot.updateFrameKinematics(q_pinocchio);
    const Eigen::Matrix3d R = robot.frameRotation(4);
    Eigen::VectorXd v_raisim(Eigen::VectorXd::Zero(18));
    // Change base velocity's reference frame from LOCAL to WORLD.
    v_raisim.head<3>().noalias()     = R * v_pinocchio.head<3>();
    // Change base angle velocity's reference frame from LOCAL to WORLD.
    v_raisim.segment<3>(3).noalias() = R * v_pinocchio.segment<3>(3);
    // LH 
    v_raisim.coeffRef(6)             = v_pinocchio.coeff(6);
    v_raisim.coeffRef(7)             = v_pinocchio.coeff(7);
    v_raisim.coeffRef(8)             = v_pinocchio.coeff(8);
    // Change leg order LH <-> RF
    v_raisim.coeffRef(9)             = v_pinocchio.coeff(12);
    v_raisim.coeffRef(10)            = v_pinocchio.coeff(13);
    v_raisim.coeffRef(11)            = v_pinocchio.coeff(14);
    // Change leg order LH <-> RF
    v_raisim.coeffRef(12)            = v_pinocchio.coeff(9);
    v_raisim.coeffRef(13)            = v_pinocchio.coeff(10);
    v_raisim.coeffRef(14)            = v_pinocchio.coeff(11);
    // RH
    v_raisim.coeffRef(15)            = v_pinocchio.coeff(15);
    v_raisim.coeffRef(16)            = v_pinocchio.coeff(16);
    v_raisim.coeffRef(17)            = v_pinocchio.coeff(17);
    return v_raisim;
  }
  else {
    return v_pinocchio;
  }
}


Eigen::VectorXd rai2pino_v(Robot& robot, const Eigen::VectorXd& q_pinocchio, 
                           const Eigen::VectorXd& v_raisim) {
  if (v_raisim.size() == 18) { // quadruped velocity
    assert(q_pinocchio.size() === 19);
    robot.updateFrameKinematics(q_pinocchio);
    const Eigen::Matrix3d R = robot.frameRotation(4);
    Eigen::VectorXd v_pinocchio(Eigen::VectorXd::Zero(18));
    // Change base velocity's reference frame from WORLD to LOCAL.
    v_pinocchio.head<3>().noalias()     = R.transpose() * v_raisim.head<3>();
    // Change base angular velocity's reference frame from WORLD to LOCAL.
    v_pinocchio.segment<3>(3).noalias() = R.transpose() * v_raisim.segment<3>(3);
    // LH 
    v_pinocchio.coeffRef(6)             = v_raisim.coeff(6);
    v_pinocchio.coeffRef(7)             = v_raisim.coeff(7);
    v_pinocchio.coeffRef(8)             = v_raisim.coeff(8);
    // Change leg order LH <-> RF
    v_pinocchio.coeffRef(9)             = v_raisim.coeff(12);
    v_pinocchio.coeffRef(10)            = v_raisim.coeff(13);
    v_pinocchio.coeffRef(11)            = v_raisim.coeff(14);
    // Change leg order LH <-> RF
    v_pinocchio.coeffRef(12)            = v_raisim.coeff(9);
    v_pinocchio.coeffRef(13)            = v_raisim.coeff(10);
    v_pinocchio.coeffRef(14)            = v_raisim.coeff(11);
    // RH
    v_pinocchio.coeffRef(15)            = v_raisim.coeff(15);
    v_pinocchio.coeffRef(16)            = v_raisim.coeff(16);
    v_pinocchio.coeffRef(17)            = v_raisim.coeff(17);
    return v_pinocchio;
  }
  else {
    return v_raisim;
  }
}


Eigen::VectorXd pino2rai_u(const Eigen::VectorXd& u_pinocchio) {
  if (u_pinocchio.size() == 12) { // quadruped torques 
    Eigen::VectorXd u_raisim(Eigen::VectorXd::Zero(18));
      // LH 
    u_raisim.coeffRef(6)  = u_pinocchio.coeff(0);
    u_raisim.coeffRef(7)  = u_pinocchio.coeff(1);
    u_raisim.coeffRef(8)  = u_pinocchio.coeff(2);
    // Change leg order LH <-> RF
    u_raisim.coeffRef(9)  = u_pinocchio.coeff(6);
    u_raisim.coeffRef(10) = u_pinocchio.coeff(7);
    u_raisim.coeffRef(11) = u_pinocchio.coeff(8);
    // Change leg order LH <-> RF
    u_raisim.coeffRef(12) = u_pinocchio.coeff(3);
    u_raisim.coeffRef(13) = u_pinocchio.coeff(4);
    u_raisim.coeffRef(14) = u_pinocchio.coeff(5);
    // RH
    u_raisim.coeffRef(15) = u_pinocchio.coeff(9);
    u_raisim.coeffRef(16) = u_pinocchio.coeff(10);
    u_raisim.coeffRef(17) = u_pinocchio.coeff(11);
    return u_raisim;
  }
  else {
    return u_pinocchio;
  }
}

} // namespace helper
} // namespace idocp