#ifndef INVDYNOCP_CONFIGURATION_SPACE_COST_HPP_
#define INVDYNOCP_CONFIGURATION_SPACE_COST_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace invdynocp {

class ConfigurationSpaceCost {
public:
  ConfigurationSpaceCost(const Robot* robot_ptr, const Eigen::VectorXd& q_weight,  
                         const Eigen::VectorXd& v_weight, 
                         const Eigen::VectorXd& a_weight,  
                         const Eigen::VectorXd& u_weight);

  ConfigurationSpaceCost(const Robot* robot_ptr, const Eigen::VectorXd& q_ref, 
                         const Eigen::VectorXd& q_weight, 
                         const Eigen::VectorXd& v_ref, 
                         const Eigen::VectorXd& v_weight, 
                         const Eigen::VectorXd& a_ref, 
                         const Eigen::VectorXd& a_weight, 
                         const Eigen::VectorXd& u_ref, 
                         const Eigen::VectorXd& u_weight);

  void lq(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& q, 
          Eigen::VectorXd& lq);

  void lv(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& v, 
          Eigen::VectorXd& lv);

  void la(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& a, 
          Eigen::VectorXd& la);

  void lu(const Robot* robot_ptr, const double dtau, const Eigen::VectorXd& u, 
          Eigen::VectorXd& lu);

  void lqq(const Robot* robot_ptr, const double dtau, Eigen::MatrixXd& lqq);

  void lvv(const Robot* robot_ptr, const double dtau, Eigen::MatrixXd& lvv);

  void laa(const Robot* robot_ptr, const double dtau, Eigen::MatrixXd& laa);

  void luu(const Robot* robot_ptr, const double dtau, Eigen::MatrixXd& luu);

  void phiq(const Robot* robot_ptr, const Eigen::VectorXd& q, 
            Eigen::VectorXd& phiq);

  void phiv(const Robot* robot_ptr, const Eigen::VectorXd& v, 
            Eigen::VectorXd& phiv);

  void phiqq(const Robot* robot_ptr, Eigen::MatrixXd& phiqq);

  void phivv(const Robot* robot_ptr, Eigen::MatrixXd& phivv);

  // Prohibits copy constructor.
  ConfigurationSpaceCost(const ConfigurationSpaceCost&) = delete;

  // Prohibits copy operator.
  ConfigurationSpaceCost& operator=(const ConfigurationSpaceCost&) = delete;

private:
  unsigned int dimq_, dimv_;
  Eigen::VectorXd q_ref_, v_ref_, a_ref_, u_ref_, q_weight_, v_weight_, 
                  a_weight_, u_weight_;

};

} // namespace invdynocp


#endif // INVDYNOCP_CONFIGURATION_SPACE_COST_HPP_