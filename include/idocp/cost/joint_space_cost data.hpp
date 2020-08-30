#ifndef IDOCP_JOINT_SPACE_COST_DATA_HPP_ 
#define IDOCP_JOINT_SPACE_COST_DATA_HPP_ 

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/cost/cost_function_component_data_base.hpp"


namespace idocp {

class JointSpaceCostData final : public CostFunctionComponentDataBase {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  JointSpaceCostData(const Robot& robot) {
    qdiff(),
    J_qdiff() {
    if (robot.has_floating_base()) {
      qdiff.resize(robot.dimv());
      qdiff.setZero();
      J_qdiff.resize(robot.dimv(), robot.dimv());
      J_qdiff.setZero();
    }
  }

  JointSpaceCostData() {
    qdiff(),
    J_qdiff() {
  }

  ~JointSpaceCostData();

  // Use defalut copy constructor.
  JointSpaceCostData(const JointSpaceCostData&) = default;

  // Use defalut copy operator.
  JointSpaceCostData& operator=(const JointSpaceCostData&) = default;

  // Use defalut move constructor.
  JointSpaceCostData(JointSpaceCostData&&) noexcept = default;

  // Use defalut move assign operator.
  JointSpaceCostData& operator=(JointSpaceCostData&&) noexcept = default;

  Eigen::VectorXd qdiff;
  Eigen::MatrixXd J_qdiff;
};

} // namespace idocp


#endif // IDOCP_JOINT_SPACE_COST_DATA_HPP_