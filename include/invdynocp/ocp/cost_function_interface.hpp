#ifndef INVDYNOCP_COST_FUNCTION_INTERFACE_HPP_
#define INVDYNOCP_COST_FUNCTION_INTERFACE_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace invdynocp {

class CostFunctionInterface {
public:
  CostFunctionInterface() {}

  virtual ~CostFunctionInterface() {}

  virtual void lq(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lq) = 0;

  virtual void lq(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::VectorXd& lq) = 0;

  virtual void lv(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lv) = 0;

  virtual void lv(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::VectorXd& lv) = 0;

  virtual void la(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& la) = 0;

  virtual void la(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::VectorXd& la) = 0;

  virtual void lu(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lu) = 0;

  virtual void lu(const Robot* robot_ptr, const double t, 
                  const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                  const Eigen::VectorXd& fext, Eigen::VectorXd& lu) = 0;

  virtual void lfext(const Robot* robot_ptr, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                     Eigen::VectorXd& lfext) = 0;

  virtual void lfext(const Robot* robot_ptr, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                     const Eigen::VectorXd& fext, Eigen::VectorXd& lfext) = 0;

  virtual void phiq(const Robot* robot_ptr, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiq) = 0;

  virtual void phiv(const Robot* robot_ptr, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiv) = 0;

  // Prohibits copy constructor.
  CostFunctionInterface(const CostFunctionInterface&) = delete;

  // Prohibits copy operator.
  CostFunctionInterface& operator=(const CostFunctionInterface&) = delete;

};

} // namespace invdynocp

#endif // INVDYNOCP_COST_FUNCTION_INTERFACE_HPP_