#ifndef IDOCP_COST_FUNCTION_INTERFACE_HPP_
#define IDOCP_COST_FUNCTION_INTERFACE_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class CostFunctionInterface {
public:
  CostFunctionInterface() {}

  virtual ~CostFunctionInterface() {}

  virtual double l(const Robot& robot, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, const Eigen::VectorXd& u) = 0;

  virtual double phi(const Robot& robot, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v) = 0;

  virtual void lq(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, Eigen::VectorXd& lq) = 0;

  virtual void lv(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, Eigen::VectorXd& lv) = 0;

  virtual void la(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, Eigen::VectorXd& la) = 0;

  virtual void lu(const Robot& robot, const double t, const double dtau,
                  const Eigen::VectorXd& u, Eigen::VectorXd& lu) = 0;

  virtual void lqq(const Robot& robot, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, Eigen::MatrixXd& lqq) = 0;

  virtual void lvv(const Robot& robot, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, Eigen::MatrixXd& lvv) = 0;

  virtual void laa(const Robot& robot, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, Eigen::MatrixXd& laa) = 0;

  virtual void luu(const Robot& robot, const double t, const double dtau,
                   const Eigen::VectorXd& u, Eigen::MatrixXd& laa) = 0;

  virtual void phiq(const Robot& robot, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiq) = 0;

  virtual void phiv(const Robot& robot, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiv) = 0;

  virtual void phiqq(const Robot& robot, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     Eigen::MatrixXd& phiqq) = 0;

  virtual void phivv(const Robot& robot, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     Eigen::MatrixXd& phivv) = 0;

  // Prohibits copy constructor.
  CostFunctionInterface(const CostFunctionInterface&) = delete;

  // Prohibits copy operator.
  CostFunctionInterface& operator=(const CostFunctionInterface&) = delete;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_INTERFACE_HPP_