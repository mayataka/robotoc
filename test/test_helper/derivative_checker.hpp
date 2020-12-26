#ifndef IDOCP_DERIVATIVE_CHECKER_HPP_
#define IDOCP_DERIVATIVE_CHECKER_HPP_

#include <memory>

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/cost/cost_function_component_base.hpp"


namespace idocp {

class DerivativeChecker {
public:
  explicit DerivativeChecker(const Robot& robot, 
                             const double finite_diff=1.0e-08, 
                             const double test_tol=0.0001);

  ~DerivativeChecker();

  void setFiniteDifference(const double finite_diff=1.0e-08);

  void setTestTolerance(const double test_tol=1.0e-08);

  bool checkFirstOrderStageCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost, 
      const ContactStatus& contact_status);

  bool checkSecondOrderStageCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost,
      const ContactStatus& contact_status);

  bool checkFirstOrderTerminalCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost);

  bool checkSecondOrderTerminalCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost);

  bool checkFirstOrderImpulseCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost, 
      const ImpulseStatus& impulse_status);

  bool checkSecondOrderImpulseCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost,
      const ImpulseStatus& impulse_status);

private:  
  Robot robot_;
  double finite_diff_, test_tol_;

};
  
} // namespace idocp 

#endif // IDOCP_DERIVATIVE_CHECKER_HPP_