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
      qdiff_3d(Eigen::Vector3d::Zero()),
      J_6d(Eigen::MatrixXd::Zero(6, robot.dimv())),
      J_3d(Eigen::MatrixXd::Zero(3, robot.dimv())) {
    if (robot.has_floating_base()) {
      qdiff.resize(robot.dimv());
      qdiff.setZero();
      J_qdiff.resize(robot.dimv(), robot.dimv());
      J_qdiff.setZero();
    }
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
  Eigen::MatrixXd J_6d, J_3d;

private:

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_DATA_HPP_ 