#ifndef IDOCP_COST_FUNCTION_DATA_HPP_
#define IDOCP_COST_FUNCTION_DATA_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class CostFunctionData
/// @brief Composed of data used to compute the cost function and its 
/// derivatives. 
///
class CostFunctionData {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW 

  ///
  /// @brief Constructor. 
  /// @param[in] robot Robot model. 
  ///
  CostFunctionData(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  CostFunctionData();

  ///
  /// @brief Destructor. 
  ///
  ~CostFunctionData();

  ///
  /// @brief Default copy constructor. 
  ///
  CostFunctionData(const CostFunctionData&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  CostFunctionData& operator=(const CostFunctionData&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  CostFunctionData(CostFunctionData&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  CostFunctionData& operator=(CostFunctionData&&) noexcept = default;

  ///
  /// @brief Vector used for computing the difference of the configurations in
  /// JointSpaceCost. 
  /// Be allocated only when Robot::has_floating_base() is true. Then the size 
  /// is Robot::dimv().
  ///
  Eigen::VectorXd qdiff;

  ///
  /// @brief Vector used for computing the difference of the position of the 
  /// end-effector in TaskSpace3DCost. 
  /// Be allocated only when CostFunction has TaskSpace3DCost. Then the size 
  /// is 3.
  ///
  Eigen::VectorXd diff_3d;

  ///
  /// @brief Vector used for computing the difference of the position and 
  /// rotation of the end-effector in TaskSpace6DCost. 
  /// Be allocated only when CostFunction has TaskSpace6DCost. Then the size 
  /// is 6.
  ///
  Eigen::VectorXd diff_6d;

  ///
  /// @brief Vector used for computing the difference of the SE3 of the 
  /// end-effector in TaskSpace6DCost. 
  /// Be allocated only when CostFunction has TaskSpace6DCost. 
  ///
  pinocchio::SE3 diff_SE3;

  ///
  /// @brief Jacobian of the difference of the configurations used in 
  /// JointSpaceCost. 
  /// Be allocated only when Robot::has_floating_base() is true. Then the size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::MatrixXd J_qdiff;

  ///
  /// @brief Jacobian used in TaskSpace3DCost and TaskSpace6DCost.
  /// Be allocated only when CostFunction has TaskSpace3DCost or 
  /// TaskSpace6DCost. Size is 6 x Robot::dimv().
  ///
  Eigen::MatrixXd J_6d;

  ///
  /// @brief Jacobian used in TaskSpace3DCost and TaskSpace6DCost.
  /// Be allocated only when CostFunction has TaskSpace3DCost or 
  /// TaskSpace6DCost. Size is 3 x Robot::dimv().
  ///
  Eigen::MatrixXd J_3d;

  ///
  /// @brief Jacobian used in TaskSpace6DCost.  Be allocated only when 
  /// CostFunction has TaskSpace6DCost. Size is 6 x 6.
  ///
  Eigen::MatrixXd J_66;

  ///
  /// @brief Jacobian used in TaskSpace6DCost. Be allocated only when 
  /// CostFunction has TaskSpace6DCost. Size is 6 x Robot::dimv().
  ///
  Eigen::MatrixXd JJ_6d;
};

} // namespace idocp

#include "idocp/cost/cost_function_data.hxx"

#endif // IDOCP_COST_FUNCTION_DATA_HPP_ 