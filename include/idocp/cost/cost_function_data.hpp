#ifndef IDOCP_COST_FUNCTION_DATA_HPP_
#define IDOCP_COST_FUNCTION_DATA_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class CostFunctionData {
public:
  CostFunctionData(const Robot& robot) 
    : qdiff(),
      J_qdiff(),
      qdiff_3d(),
      J_6d(),
      J_3d() {
    if (robot.has_floating_base()) {
      qdiff.resize(robot.dimv());
      qdiff.setZero();
      J_qdiff.resize(robot.dimv(), robot.dimv());
      J_qdiff.setZero();
    }
    J_6d.setZero();
    J_3d.setZero();
    qdiff_3d.setZero();
  }

  CostFunctionData() 
    : qdiff(),
      J_qdiff(),
      qdiff_3d(),
      J_6d(),
      J_3d() {
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

  Eigen::VectorXd qdiff;
  Eigen::MatrixXd J_qdiff;
  Eigen::Vector3d qdiff_3d;
  Eigen::Matrix<double, 6, Eigen::Dynamic> J_6d;
  Eigen::Matrix<double, 3, Eigen::Dynamic> J_3d;

private:

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_DATA_HPP_ 