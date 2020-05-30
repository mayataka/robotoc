#ifndef INVDYNOCP_CONSTRAINTS_INTERFACE_HPP_
#define INVDYNOCP_CONSTRAINTS_INTERFACE_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace invdynocp {

class ConstraintsInterface {
public:
  ConstraintsInterface() {}

  virtual ~ConstraintsInterface() {}

  virtual void C(const Robot* robot_ptr, const double t, 
                 const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                 const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                 Eigen::VectorXd& C) = 0;

  virtual void Cq(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Cq) = 0;

  virtual void Cq(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Cq) = 0;

  virtual void Cv(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Cv) = 0;

  virtual void Cv(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Cv) = 0;

  virtual void Ca(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Ca) = 0;

  virtual void Ca(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Ca) = 0;

  virtual void Cu(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::MatrixXd& Cu) = 0;

  virtual void Cu(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::MatrixXd& Cu) = 0;

  virtual void Cfext(const Robot* robot_ptr, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                     Eigen::MatrixXd& Cfext) = 0;

  virtual void Cfext(const Robot* robot_ptr, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                     const Eigen::VectorXd& fext, Eigen::MatrixXd& Cfext) = 0;

  virtual unsigned int dimc() const = 0;

  // Prohibits copy constructor.
  ConstraintsInterface(const ConstraintsInterface&) = delete;

  // Prohibits copy operator.
  ConstraintsInterface& operator=(const ConstraintsInterface&) = delete;

};

} // namespace invdynocp

#endif // INVDYNOCP_CONSTRAINTS_INTERFACE_HPP_