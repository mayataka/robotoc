#ifndef IDOCP_COST_FUNCTION_INTERFACE_HPP_
#define IDOCP_COST_FUNCTION_INTERFACE_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class CostFunctionInterface {
public:
  CostFunctionInterface() {}

  virtual ~CostFunctionInterface() {}

  // Use default copy constructor.
  CostFunctionInterface(const CostFunctionInterface&) = default;

  // Use default copy coperator.
  CostFunctionInterface& operator=(const CostFunctionInterface&) = default;

  virtual void setConfigurationJacobian(const Robot& robot, 
                                        const Eigen::VectorXd& q) = 0;

  virtual double l(const double t, const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   const Eigen::VectorXd& u) = 0;

  virtual double phi(const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v) = 0;

  virtual void lq(const double t, const double dtau, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  Eigen::VectorXd& lq) = 0;

  virtual void lv(const double t, const double dtau, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  Eigen::VectorXd& lv) = 0;

  virtual void la(const double t, const double dtau, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  Eigen::VectorXd& la) = 0;

  virtual void lu(const double t, const double dtau, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& lu) = 0;

  virtual void lf(const double t, const double dtau, const Eigen::VectorXd& f, 
                  Eigen::VectorXd& lf) = 0;

  virtual void lqq(const double t, const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lqq) = 0;

  virtual void lvv(const double t, const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lvv) = 0;

  virtual void laa(const double t, const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& laa) = 0;

  virtual void luu(const double t, const double dtau, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& luu) = 0;

  virtual void lff(const double t, const double dtau, const Eigen::VectorXd& f, 
                   Eigen::MatrixXd& lff) = 0;

  virtual void augment_lqq(const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, Eigen::MatrixXd& lqq) = 0;

  virtual void augment_lvv(const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, Eigen::MatrixXd& lvv) = 0;

  virtual void augment_laa(const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, Eigen::MatrixXd& laa) = 0;

  virtual void augment_luu(const double t, const double dtau, 
                           const Eigen::VectorXd& u, Eigen::MatrixXd& luu) = 0;

  virtual void augment_lff(const double t, const double dtau, 
                           const Eigen::VectorXd& f, Eigen::MatrixXd& lff) = 0;

  virtual void phiq(const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::VectorXd& phiq) = 0;

  virtual void phiv(const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::VectorXd& phiv) = 0;

  virtual void phiqq(const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v, Eigen::MatrixXd& phiqq) = 0;

  virtual void phivv(const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v, Eigen::MatrixXd& phivv) = 0;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_INTERFACE_HPP_