#include "constraints/configuration_space_constraints.hpp"


namespace invdynocp {


ConfigurationSpaceConstraints::ConfigurationSpaceConstraints(
    const Robot* robot_ptr) {
}


ConfigurationSpaceConstraints::ConfigurationSpaceConstraints(
    const Robot* robot_ptr, const Eigen::VectorXd& q_max, 
    const Eigen::VectorXd& q_min, const Eigen::VectorXd& v_max, 
    const Eigen::VectorXd& v_min, const Eigen::VectorXd& u_max, 
    const Eigen::VectorXd& u_min)
  : q_max_(q_max),
    q_min_(q_min),
    v_max_(v_max),
    v_min_(v_min),
    u_max_(u_max),
    u_min_(u_min) {
}


void ConfigurationSpaceConstraints::C(const Robot* robot_ptr, const double dtau, 
                                      const Eigen::VectorXd& q, 
                                      const Eigen::VectorXd& v, 
                                      const Eigen::VectorXd& u, 
                                      Eigen::VectorXd& C) {
}


void ConfigurationSpaceConstraints::Cq(const Robot* robot_ptr, 
                                       const double dtau, 
                                       const Eigen::VectorXd& q, 
                                       Eigen::MatrixXd& Cq) {
}


void ConfigurationSpaceConstraints::Cv(const Robot* robot_ptr, 
                                       const double dtau, 
                                       const Eigen::VectorXd& v, 
                                       Eigen::MatrixXd& Cv) {
}


void ConfigurationSpaceConstraints::Cu(const Robot* robot_ptr, 
                                       const double dtau, 
                                       const Eigen::VectorXd& u, 
                                       Eigen::MatrixXd& Cu) {
}

}