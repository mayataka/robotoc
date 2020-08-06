#ifndef IDOCP_COST_FUNCTION_INTERFACE_HPP_
#define IDOCP_COST_FUNCTION_INTERFACE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class CostFunctionInterface {
public:
  CostFunctionInterface() {}

  virtual ~CostFunctionInterface() {}

  // Use default copy constructor.
  CostFunctionInterface(const CostFunctionInterface&) = default;

  // Use default copy coperator.
  CostFunctionInterface& operator=(const CostFunctionInterface&) = default;

  // Use default move constructor.
  CostFunctionInterface(CostFunctionInterface&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionInterface& operator=(CostFunctionInterface&&) noexcept = default;

  virtual double l(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   const Eigen::VectorXd& u, const Eigen::VectorXd& f) const = 0;

  virtual double phi(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::VectorXd& q, 
                     const Eigen::VectorXd& v) const = 0;

  virtual void lq(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  Eigen::VectorXd& lq) const = 0;

  virtual void lv(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  Eigen::VectorXd& lv) const = 0;

  virtual void la(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                  Eigen::VectorXd& la) const = 0;

  virtual void lu(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& la) const = 0;

  virtual void lf(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::VectorXd& u, 
                  Eigen::VectorXd& la) const = 0;

  virtual void lqq(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lqq) const = 0;

  virtual void lvv(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& lvv) const = 0;

  virtual void laa(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& q, 
                   const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                   Eigen::MatrixXd& laa) const = 0;

  virtual void luu(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& u, 
                   Eigen::MatrixXd& luu) const = 0;

  virtual void lff(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, const Eigen::VectorXd& f, 
                   Eigen::MatrixXd& lff) const = 0;

  virtual void augment_lqq(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, 
                           Eigen::MatrixXd& lqq) const = 0;

  virtual void augment_lvv(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, 
                           Eigen::MatrixXd& lvv) const = 0;

  virtual void augment_laa(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                           const Eigen::VectorXd& a, 
                           Eigen::MatrixXd& laa) const = 0;

  virtual void augment_luu(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::VectorXd& u, 
                           Eigen::MatrixXd& luu) const = 0;

  virtual void augment_lff(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::VectorXd& f, 
                           Eigen::MatrixXd& lff) const = 0;

  virtual void phiq(const Robot& robot, CostFunctionData& data, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiq) const = 0;

  virtual void phiv(const Robot& robot, CostFunctionData& data, const double t, 
                    const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                    Eigen::VectorXd& phiv) const = 0;

  virtual void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     Eigen::MatrixXd& phiqq) const = 0;

  virtual void phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                     Eigen::MatrixXd& phivv) const = 0;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_INTERFACE_HPP_