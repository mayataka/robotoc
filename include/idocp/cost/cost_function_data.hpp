#ifndef IDOCP_COST_FUNCTION_DATA_HPP_
#define IDOCP_COST_FUNCTION_DATA_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class CostFunctionData {
public:
  CostFunctionData(const Robot& robot);

  CostFunctionData();

  ~CostFunctionData();

  // Use default copy constructor.
  CostFunctionData(const CostFunctionData&) = default;

  // Use default copy coperator.
  CostFunctionData& operator=(const CostFunctionData&) = default;

  // Use default move constructor.
  CostFunctionData(CostFunctionData&&) noexcept = default;

  // Use default move assign coperator.
  CostFunctionData& operator=(CostFunctionData&&) noexcept = default;

  Eigen::VectorXd qdiff, diff_3d, diff_6d;
  pinocchio::SE3 diff_SE3;
  Eigen::MatrixXd J_qdiff;
  Eigen::MatrixXd J_6d, J_3d;
  Eigen::MatrixXd J_66, JJ_6d;
};

} // namespace idocp

#include "idocp/cost/cost_function_data.hxx"

#endif // IDOCP_COST_FUNCTION_DATA_HPP_ 