#ifndef IDOCP_CONFIGURATION_SPACE_CONSTRAINTS_HPP_
#define IDOCP_CONFIGURATION_SPACE_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class ConfigurationSpaceConstraints {
public:
  ConfigurationSpaceConstraints(const Robot* robot_ptr);

  ConfigurationSpaceConstraints(const Robot* robot_ptr, 
                                const Eigen::VectorXd& q_max, 
                                const Eigen::VectorXd& q_min, 
                                const Eigen::VectorXd& v_max, 
                                const Eigen::VectorXd& v_min, 
                                const Eigen::VectorXd& u_max,
                                const Eigen::VectorXd& u_min);

  void C(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& q, 
         const Eigen::VectorXd& v, const Eigen::VectorXd& u, 
         Eigen::VectorXd& C);

  void Cq(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& q, 
          Eigen::MatrixXd& Cq);

  void Cv(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& v, 
          Eigen::MatrixXd& Cv);

  void Cu(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& u, 
          Eigen::MatrixXd& Cu);

  // Prohibits copy constructor.
  ConfigurationSpaceConstraints(const ConfigurationSpaceConstraints&) = delete;

  // Prohibits copy operator.
  ConfigurationSpaceConstraints& operator=(const ConfigurationSpaceConstraints&) 
      = delete;

private:
  Eigen::VectorXd q_max_, q_min_, v_max_, v_min_, u_max_, u_min_;

};

} // namespace idocp

#endif // IDOCP_CONFIGURATION_SPACE_CONSTRAINTS_HPP_