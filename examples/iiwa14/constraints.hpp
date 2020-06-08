#ifndef INVDYNOCP_IIWA14_CONSTRAINTS_HPP_
#define INVDYNOCP_IIWA14_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/constraints_interface.hpp"
#include "constraints/configuration_space_constraints.hpp"


namespace invdynocp {
namespace iiwa14 {

class Constraints final : public ConstraintsInterface {
public:
  Constraints(const Robot* robot_ptr);

  void C(const Robot* robot_ptr, const double t, const double dtau, 
         const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
         const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
         Eigen::VectorXd& C) override;

  void Cq(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          Eigen::MatrixXd& Cq) override;

  void Cq(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          const Eigen::VectorXd& fext, Eigen::MatrixXd& Cq) override;

  void Cv(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          Eigen::MatrixXd& Cv) override;

  void Cv(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          const Eigen::VectorXd& fext, Eigen::MatrixXd& Cv) override;

  void Ca(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          Eigen::MatrixXd& Ca) override;

  void Ca(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          const Eigen::VectorXd& fext, Eigen::MatrixXd& Ca) override;

  void Cu(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          Eigen::MatrixXd& Cu) override;

  void Cu(const Robot* robot_ptr, const double t, const double dtau,
          const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
          const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
          const Eigen::VectorXd& fext, Eigen::MatrixXd& Cu) override;

  void Cfext(const Robot* robot_ptr, const double t, const double dtau,
             const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
             const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
             const Eigen::VectorXd& fext, Eigen::MatrixXd& Cfext) override;

  unsigned int dimc() const override;


private:
  ConfigurationSpaceConstraints configuration_space_constraints_;  

};

} // namespace iiwa14
} // namespace invdynocp

#endif // INVDYNOCP_IIWA14_CONSTRAINTS_HPP_