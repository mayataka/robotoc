#ifndef IDOCP_JOINT_TORQUE_UPPER_LIMITS_HPP_
#define IDOCP_JOINT_TORQUE_UPPER_LIMITS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class JointTorqueUpperLimits {
public:
  JointTorqueUpperLimits(const Robot& robot);

  // Use default copy constructor.
  JointTorqueUpperLimits(const JointTorqueUpperLimits&) = default;

  // Use default copy operator.
  JointTorqueUpperLimits& operator=(const JointTorqueUpperLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double barrier, 
                       const Eigen::VectorXd& u);

  void linearizeConstraint(const Robot& robot, const double barrier, 
                           const Eigen::VectorXd& u);

  void updateSlackAndDual(const Robot& robot, 
                          const Eigen::MatrixXd& du_dq,
                          const Eigen::MatrixXd& du_dv,
                          const Eigen::MatrixXd& du_da,
                          const Eigen::VectorXd& dq,
                          const Eigen::VectorXd& dv,
                          const Eigen::VectorXd& du); 

  void condenseSlackAndDual(const Robot& robot, const Eigen::MatrixXd& du_dq,
                            const Eigen::MatrixXd& du_dv,
                            const Eigen::MatrixXd& du_da, Eigen::MatrixXd& Cqq, 
                            Eigen::MatrixXd& Cqv, Eigen::MatrixXd& Cqa, 
                            Eigen::MatrixXd& Cvv, Eigen::MatrixXd& Cva, 
                            Eigen::MatrixXd& Caa, Eigen::VectorXd& Cq,
                            Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void augmentDualResidual(const Robot& robot, const Eigen::MatrixXd& du_dq,
                           const Eigen::MatrixXd& du_dv, 
                           const Eigen::MatrixXd& du_da, Eigen::VectorXd& Cq,
                           Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  double squaredConstraintsResidualNrom(const Robot& robot, 
                                        const Eigen::VectorXd& u);

private:
  unsigned int dimq_, dimv_, dimc_;
  Eigen::VectorXd umax_, slack_, dual_, residual_, FB_residual_, dslack_, 
                  ddual_, newton_residual_;
  Eigen::MatrixXd partial_du_;
};

} // namespace idocp


#endif // IDOCP_JOINT_TORQUE_UPPER_LIMITS_HPP_