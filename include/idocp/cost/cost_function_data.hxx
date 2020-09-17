#ifndef IDOCP_COST_FUNCTION_DATA_HXX_
#define IDOCP_COST_FUNCTION_DATA_HXX_

#include "idocp/cost/cost_function_data.hpp"

namespace idocp {

inline CostFunctionData::CostFunctionData(const Robot& robot) 
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


inline CostFunctionData::CostFunctionData() 
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


inline CostFunctionData::~CostFunctionData() {
}

} // namespace idocp


#endif // IDOCP_COST_FUNCTION_DATA_HXX_ 