#ifndef IDOCP_COST_FUNCTION_DATA_HPP_
#define IDOCP_COST_FUNCTION_DATA_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class CostFunctionData {
public:
  CostFunctionData(const Robot& robot) 
    : qdiff(),
      diff_3d(Eigen::VectorXd::Zero(3)),
      diff_6d(Eigen::VectorXd::Zero(6)),
      diff_SE3(pinocchio::SE3(Eigen::Matrix3d::Identity(), 
                              Eigen::Vector3d::Zero())),
      J_qdiff(),
      J_6d(Eigen::MatrixXd::Zero(6, robot.dimv())),
      J_3d(Eigen::MatrixXd::Zero(3, robot.dimv())),
      J_66(Eigen::MatrixXd::Zero(6, 6)),
      JJ_6d(Eigen::MatrixXd::Zero(6, robot.dimv())) {
    if (robot.has_floating_base()) {
      qdiff.resize(robot.dimv());
      qdiff.setZero();
      J_qdiff.resize(robot.dimv(), robot.dimv());
      J_qdiff.setZero();
    }
  }

  CostFunctionData() 
    : qdiff(),
      diff_3d(),
      diff_6d(),
      diff_SE3(),
      J_qdiff(),
      J_6d(),
      J_3d(),
      J_66(),
      JJ_6d() {
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

  Eigen::VectorXd qdiff, diff_3d, diff_6d;
  pinocchio::SE3 diff_SE3;
  Eigen::MatrixXd J_qdiff;
  Eigen::MatrixXd J_6d, J_3d;
  Eigen::MatrixXd J_66, JJ_6d;

private:

};

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_DATA_HPP_ 