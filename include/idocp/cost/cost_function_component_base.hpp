#ifndef IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_
#define IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_data.hpp"


namespace idocp {

class CostFunctionComponentBase {
public:
  CostFunctionComponentBase() {}

  virtual ~CostFunctionComponentBase() {}

  // Use default copy constructor.
  CostFunctionComponentBase(const CostFunctionComponentBase&) = default;

  // Use default copy coperator.
  CostFunctionComponentBase& operator=(const CostFunctionComponentBase&) 
      = default;

  // Use default move constructor.
  CostFunctionComponentBase(CostFunctionComponentBase&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionComponentBase& operator=(CostFunctionComponentBase&&) noexcept 
      = default;

  virtual double l(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   const Eigen::Ref<const Eigen::VectorXd>& f, 
                   const Eigen::Ref<const Eigen::VectorXd>& u) const = 0;

  virtual double phi(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v) const = 0;

  virtual void lq(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
                  const Eigen::Ref<const Eigen::VectorXd>& v, 
                  const Eigen::Ref<const Eigen::VectorXd>& a, 
                  Eigen::Ref<Eigen::VectorXd> lq) const = 0;

  virtual void lv(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
                  const Eigen::Ref<const Eigen::VectorXd>& v, 
                  const Eigen::Ref<const Eigen::VectorXd>& a, 
                  Eigen::Ref<Eigen::VectorXd> lv) const = 0;

  virtual void la(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::Ref<const Eigen::VectorXd>& q, 
                  const Eigen::Ref<const Eigen::VectorXd>& v, 
                  const Eigen::Ref<const Eigen::VectorXd>& a, 
                  Eigen::Ref<Eigen::VectorXd> la) const = 0;

  virtual void lf(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::Ref<const Eigen::VectorXd>& f, 
                  Eigen::Ref<Eigen::VectorXd> lf) const = 0;

  virtual void lu(const Robot& robot, CostFunctionData& data, const double t, 
                  const double dtau, const Eigen::Ref<const Eigen::VectorXd>& u, 
                  Eigen::Ref<Eigen::VectorXd> lu) const = 0;

  virtual void lqq(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lqq) const = 0;

  virtual void lvv(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> lvv) const = 0;

  virtual void laa(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& q, 
                   const Eigen::Ref<const Eigen::VectorXd>& v, 
                   const Eigen::Ref<const Eigen::VectorXd>& a, 
                   Eigen::Ref<Eigen::MatrixXd> laa) const = 0;

  virtual void lff(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& f, 
                   Eigen::Ref<Eigen::MatrixXd> lff) const = 0;

  virtual void luu(const Robot& robot, CostFunctionData& data, const double t, 
                   const double dtau, 
                   const Eigen::Ref<const Eigen::VectorXd>& u, 
                   Eigen::Ref<Eigen::MatrixXd> luu) const = 0;

  virtual void augment_lqq(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::Ref<const Eigen::VectorXd>& q, 
                           const Eigen::Ref<const Eigen::VectorXd>& v, 
                           const Eigen::Ref<const Eigen::VectorXd>& a, 
                           Eigen::Ref<Eigen::MatrixXd> lqq) const = 0;

  virtual void augment_lvv(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::Ref<const Eigen::VectorXd>& q, 
                           const Eigen::Ref<const Eigen::VectorXd>& v, 
                           const Eigen::Ref<const Eigen::VectorXd>& a, 
                           Eigen::Ref<Eigen::MatrixXd> lvv) const = 0;


  virtual void augment_laa(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::Ref<const Eigen::VectorXd>& q, 
                           const Eigen::Ref<const Eigen::VectorXd>& v, 
                           const Eigen::Ref<const Eigen::VectorXd>& a, 
                           Eigen::Ref<Eigen::MatrixXd> laa) const = 0;

  virtual void augment_lff(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::Ref<const Eigen::VectorXd>& f, 
                           Eigen::Ref<Eigen::MatrixXd> lff) const = 0;

  virtual void augment_luu(const Robot& robot, CostFunctionData& data, 
                           const double t, const double dtau, 
                           const Eigen::Ref<const Eigen::VectorXd>& u, 
                           Eigen::Ref<Eigen::MatrixXd> luu) const = 0;

  virtual void phiq(const Robot& robot, CostFunctionData& data, const double t, 
                    const Eigen::Ref<const Eigen::VectorXd>& q, 
                    const Eigen::Ref<const Eigen::VectorXd>& v, 
                    Eigen::Ref<Eigen::VectorXd> phiq) const = 0;

  virtual void phiv(const Robot& robot, CostFunctionData& data, const double t, 
                    const Eigen::Ref<const Eigen::VectorXd>& q, 
                    const Eigen::Ref<const Eigen::VectorXd>& v, 
                    Eigen::Ref<Eigen::VectorXd> phiv) const = 0;

  virtual void phiqq(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phiqq) const = 0;

  virtual void phivv(const Robot& robot, CostFunctionData& data, const double t, 
                     const Eigen::Ref<const Eigen::VectorXd>& q, 
                     const Eigen::Ref<const Eigen::VectorXd>& v, 
                     Eigen::Ref<Eigen::MatrixXd> phivv) const = 0;

  virtual void augment_phiqq(const Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const Eigen::Ref<const Eigen::VectorXd>& q, 
                             const Eigen::Ref<const Eigen::VectorXd>& v, 
                             Eigen::Ref<Eigen::MatrixXd> phiqq) const = 0;

  virtual void augment_phivv(const Robot& robot, CostFunctionData& data, 
                             const double t, 
                             const Eigen::Ref<const Eigen::VectorXd>& q, 
                             const Eigen::Ref<const Eigen::VectorXd>& v, 
                             Eigen::Ref<Eigen::MatrixXd> phivv) const = 0;

};

} // namespace idocp

#endif // IDOCP_COST_FUNCTION_COMPONENT_BASE_HPP_