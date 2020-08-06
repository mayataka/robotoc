#ifndef IDOCP_COST_FUNCTION_DATA_HPP_
#define IDOCP_COST_FUNCTION_DATA_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class CostFunctionData {
public:
  CostFunctionData(const Robot& robot) 
    : lq_configuration(),
      configuration_jacobian() {
    if (robot.has_floating_base()) {
      lq_configuration.resize(robot.dimq());
      lq_configuration.setZero();
      configuration_jacobian.resize(robot.dimq(), robot.dimv());
      configuration_jacobian.setZero();
    }
  }

  CostFunctionData() 
    : lq_configuration(),
      configuration_jacobian() {
  }

  ~CostFunctionData() {
  }

  // Use default copy constructor.
  CostFunctionData(const CostFunctionData&) = default;

  // Use default copy coperator.
  CostFunctionData& operator=(const CostFunctionData&) = default;

  // Use default move constructor.
  CostFunctionData(CostFunctionData&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionData& operator=(CostFunctionData&&) noexcept = default;

  Eigen::VectorXd lq_configuration;
  Eigen::MatrixXd configuration_jacobian;

private:

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_DATA_HPP_ 