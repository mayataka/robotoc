#ifndef IDOCP_CONSTRAINTS_INTERFACE_HPP_
#define IDOCP_CONSTRAINTS_INTERFACE_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class ConstraintsInterface {
public:
  ConstraintsInterface() {}

  virtual ~ConstraintsInterface() {}

  virtual void C(const Robot& robot, const double t, const double dtau, 
                 const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                 const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                 Eigen::VectorXd& C) = 0;

  virtual void Cq(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Cq) = 0;

  virtual void Cq(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Cq) = 0;

  virtual void Cv(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Cv) = 0;

  virtual void Cv(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Cv) = 0;

  virtual void Ca(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Ca) = 0;

  virtual void Ca(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Ca) = 0;

  virtual void Cu(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Cu) = 0;

  virtual void Cu(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Cu) = 0;

  virtual void Cfext(const Robot& robot, const double t, const double dtau,
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                     const Eigen::VectorXd& fext, Eigen::MatrixXd& Cfext) = 0;

  virtual unsigned int dimc() const = 0;

  // Prohibits copy constructor.
  ConstraintsInterface(const ConstraintsInterface&) = delete;

  // Prohibits copy operator.
  ConstraintsInterface& operator=(const ConstraintsInterface&) = delete;

};

} // namespace idocp

#endif // IDOCP_CONSTRAINTS_INTERFACE_HPP_