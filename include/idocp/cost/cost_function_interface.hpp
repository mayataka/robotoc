#ifndef IDOCP_COST_FUNCTION_INTERFACE_HPP_
#define IDOCP_COST_FUNCTION_INTERFACE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class CostFunctionInterface {
public:
  CostFunctionInterface() {}

  virtual ~CostFunctionInterface() {}

  // Use default copy constructor.
  CostFunctionInterface(const CostFunctionInterface&) = default;

  // Use default copy coperator.
  CostFunctionInterface& operator=(const CostFunctionInterface&) = default;

  virtual double l(const Robot& robot, const double t, const double dtau, 
                   const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                   const Eigen::VectorXd& a, const Eigen::VectorXd& u, 
                   const Eigen::VectorXd& f) = 0;

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
                  const Eigen::VectorXd& u, Eigen::VectorXd& la) = 0;

  virtual void lf(const Robot& robot, const double t, const double dtau, 
                  const Eigen::VectorXd& u, Eigen::VectorXd& la) = 0;

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
                   const Eigen::VectorXd& u, Eigen::MatrixXd& luu) = 0;

  virtual void lff(const Robot& robot, const double t, const double dtau, 
                   const Eigen::VectorXd& f, Eigen::MatrixXd& lff) = 0;

  virtual void augment_lqq(const Robot& robot, const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, Eigen::MatrixXd& lqq) = 0;

  virtual void augment_lvv(const Robot& robot, const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, Eigen::MatrixXd& lvv) = 0;

  virtual void augment_laa(const Robot& robot, const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, Eigen::MatrixXd& laa) = 0;

  virtual void augment_luu(const Robot& robot, const double t, const double dtau, 
                           const Eigen::VectorXd& u, Eigen::MatrixXd& luu) = 0;

  virtual void augment_lff(const Robot& robot, const double t, const double dtau, 
                           const Eigen::VectorXd& f, Eigen::MatrixXd& lff) = 0;

  virtual void phiq(const Robot& robot, const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::VectorXd& phiq) = 0;

  virtual void phiv(const Robot& robot, const double t, const Eigen::VectorXd& q, 
                    const Eigen::VectorXd& v, Eigen::VectorXd& phiv) = 0;

  virtual void phiqq(const Robot& robot, const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v, Eigen::MatrixXd& phiqq) = 0;

  virtual void phivv(const Robot& robot, const double t, const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v, Eigen::MatrixXd& phivv) = 0;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_INTERFACE_HPP_