#ifndef ROBOTOC_DERIVATIVE_CHECKER_HPP_
#define ROBOTOC_DERIVATIVE_CHECKER_HPP_

#include <memory>

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/cost/cost_function_component_base.hpp"


namespace robotoc {

class DerivativeChecker {
public:
  explicit DerivativeChecker(const Robot& robot, 
                             const double finite_diff=1.0e-08, 
                             const double test_tol=1.0e-03);

  ~DerivativeChecker();

  void setFiniteDifference(const double finite_diff=1.0e-08);

  void setTestTolerance(const double test_tol=1.0e-03);

  bool checkFirstOrderStageCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost);

  bool checkSecondOrderStageCostDerivatives(
      const std::shared_ptr<CostFunctionComponentBase>& cost);

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
      const std::shared_ptr<CostFunctionComponentBase>& cost);

  bool checkSecondOrderImpulseCostDerivatives(
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
  
} // namespace robotoc 

#endif // ROBOTOC_DERIVATIVE_CHECKER_HPP_