#ifndef IDOCP_COST_FUNCTION_DATA_HPP_
#define IDOCP_COST_FUNCTION_DATA_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class CostFunctionData {
public:
  CostFunctionData(const Robot& robot) 
    : q_diff(),
      Jq_diff() {
    if (robot.has_floating_base()) {
      q_diff.resize(robot.dimv());
      q_diff.setZero();
      Jq_diff.resize(robot.dimv(), robot.dimv());
      Jq_diff.setZero();
    }
  }

  CostFunctionData() 
    : q_diff(),
      Jq_diff() {
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

  Eigen::VectorXd q_diff;
  Eigen::MatrixXd Jq_diff;

private:

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_DATA_HPP_ 