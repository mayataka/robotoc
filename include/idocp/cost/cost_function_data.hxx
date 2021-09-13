#ifndef IDOCP_COST_FUNCTION_DATA_HXX_
#define IDOCP_COST_FUNCTION_DATA_HXX_

#include "idocp/cost/cost_function_data.hpp"

namespace idocp {

inline CostFunctionData::CostFunctionData(const Robot& robot) 
  : qdiff(Eigen::VectorXd::Zero(robot.dimv())),
    q_ref(Eigen::VectorXd::Zero(robot.dimq())),
    q_3d_ref(Eigen::VectorXd::Zero(3)),
    diff_3d(Eigen::VectorXd::Zero(3)),
    diff_6d(Eigen::VectorXd::Zero(6)),
    SE3_ref(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    SE3_ref_inv(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    diff_SE3(SE3(Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero())),
    J_qdiff(),
    J_6d(Eigen::MatrixXd::Zero(6, robot.dimv())),
    J_3d(Eigen::MatrixXd::Zero(3, robot.dimv())),
    J_66(Eigen::MatrixXd::Zero(6, 6)),
    JJ_6d(Eigen::MatrixXd::Zero(6, robot.dimv())) {
  if (robot.hasFloatingBase()) {
    qdiff.resize(robot.dimv());
    qdiff.setZero();
    J_qdiff.resize(robot.dimv(), robot.dimv());
    J_qdiff.setZero();
  }
}


inline CostFunctionData::CostFunctionData() 
  : qdiff(),
    q_ref(),
    q_3d_ref(),
    diff_3d(),
    diff_6d(),
    SE3_ref(),
    SE3_ref_inv(),
    diff_SE3(),
    J_qdiff(),
    J_6d(),
    J_3d(),
    J_66(),
    JJ_6d() {
}


inline CostFunctionData::~CostFunctionData() {
}

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_DATA_HXX_ 