#ifndef IDOCP_CONSTRAINTS_INTERFACE_HPP_
#define IDOCP_CONSTRAINTS_INTERFACE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class ConstraintsInterface {
public:
  ConstraintsInterface() {}

  virtual ~ConstraintsInterface() {}

  // Use default copy constructor.
  ConstraintsInterface(const ConstraintsInterface&) = default;

  // Use default copy coperator.
  ConstraintsInterface& operator=(const ConstraintsInterface&) = default;

  // Use default move constructor.
  ConstraintsInterface(ConstraintsInterface&&) noexcept = default;

  // Use default move assign coperator.
  ConstraintsInterface& operator=(ConstraintsInterface&&) noexcept = default;

  virtual bool isFeasible(const Robot& robot, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                          const Eigen::VectorXd& u) = 0;

  virtual void setSlackAndDual(const Robot& robot, const double dtau, 
                               const Eigen::VectorXd& q, 
                               const Eigen::VectorXd& v, 
                               const Eigen::VectorXd& a, 
                               const Eigen::VectorXd& u) = 0;

  virtual void augmentDualResidual(const Robot& robot, const double dtau, 
                                   Eigen::VectorXd& Cu) = 0;

  virtual void augmentDualResidual(const Robot& robot, const double dtau, 
                                   Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                                   Eigen::VectorXd& Ca) = 0;

  virtual void condenseSlackAndDual(const Robot& robot, const double dtau, 
                                    const Eigen::VectorXd& q, 
                                    const Eigen::VectorXd& v, 
                                    const Eigen::VectorXd& a, 
                                    Eigen::MatrixXd& Cqq, Eigen::MatrixXd& Cqv, 
                                    Eigen::MatrixXd& Cqa,  Eigen::MatrixXd& Cvv, 
                                    Eigen::MatrixXd& Cva,  Eigen::MatrixXd& Caa,  
                                    Eigen::VectorXd& Cq, Eigen::VectorXd& Cv, 
                                    Eigen::VectorXd& Ca) = 0;

  virtual void condenseSlackAndDual(const Robot& robot, const double dtau, 
                                    const Eigen::VectorXd& u, 
                                    Eigen::MatrixXd& Cuu, 
                                    Eigen::VectorXd& Cu) = 0;

  virtual void computeSlackAndDualDirection(const Robot& robot, 
                                            const double dtau,
                                            const Eigen::VectorXd& dq,
                                            const Eigen::VectorXd& dv,
                                            const Eigen::VectorXd& da,
                                            const Eigen::VectorXd& du) = 0;

  virtual double maxSlackStepSize() = 0;

  virtual double maxDualStepSize() = 0;

  virtual void updateSlack(const double step_size) = 0;

  virtual void updateDual(const double step_size) = 0;

  virtual double costSlackBarrier() = 0;

  virtual double costSlackBarrier(const double step_size) = 0;

  virtual double residualL1Nrom(const Robot& robot, const double dtau,
                                const Eigen::VectorXd& q, 
                                const Eigen::VectorXd& v, 
                                const Eigen::VectorXd& a, 
                                const Eigen::VectorXd& u) = 0;

  virtual double residualSquaredNrom(const Robot& robot, const double dtau,
                                     const Eigen::VectorXd& q, 
                                     const Eigen::VectorXd& v, 
                                     const Eigen::VectorXd& a, 
                                     const Eigen::VectorXd& u) = 0;

};

} // namespace idocp

#endif // IDOCP_CONSTRAINTS_INTERFACE_HPP_