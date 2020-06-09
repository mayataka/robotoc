#ifndef IDOCP_COST_FUNCTION_INTERFACE_HPP_
#define IDOCP_COST_FUNCTION_INTERFACE_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class CostFunctionInterface {
public:
  CostFunctionInterface() {}

  virtual ~CostFunctionInterface() {}

  virtual void lq(const Robot* robot_ptr, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lq) = 0;

  virtual void lv(const Robot* robot_ptr, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lv) = 0;

  virtual void la(const Robot* robot_ptr, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& la) = 0;

  virtual void lu(const Robot* robot_ptr, const double t, const double dtau,
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lu) = 0;

  virtual void lqq(const Robot* robot_ptr, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& lqq) = 0;

  virtual void lvv(const Robot* robot_ptr, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& lvv) = 0;

  virtual void laa(const Robot* robot_ptr, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& laa) = 0;

  virtual void luu(const Robot* robot_ptr, const double t, const double dtau,
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& laa) = 0;

  virtual void phiq(const Robot* robot_ptr, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiq) = 0;

  virtual void phiv(const Robot* robot_ptr, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiv) = 0;

  virtual void phiqq(const Robot* robot_ptr, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     Eigen::MatrixXd& phiqq) = 0;

  virtual void phivv(const Robot* robot_ptr, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     Eigen::MatrixXd& phivv) = 0;

  // Prohibits copy constructor.
  CostFunctionInterface(const CostFunctionInterface&) = delete;

  // Prohibits copy operator.
  CostFunctionInterface& operator=(const CostFunctionInterface&) = delete;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_INTERFACE_HPP_