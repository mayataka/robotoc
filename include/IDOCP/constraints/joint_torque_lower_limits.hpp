#ifndef IDOCP_JOINT_TORQUE_LOWER_LIMITS_HPP_
#define IDOCP_JOINT_TORQUE_LOWER_LIMITS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class JointTorqueLowerLimits {
public:
  JointTorqueLowerLimits(const Robot& robot, const double barrier);

  // Use default copy constructor.
  JointTorqueLowerLimits(const JointTorqueLowerLimits&) = default;

  // Use default copy operator.
  JointTorqueLowerLimits& operator=(const JointTorqueLowerLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double dtau, 
                       const Eigen::VectorXd& u);

  void linearizeConstraint(const Robot& robot, const double dtau, 
                           const Eigen::VectorXd& u);

  void updateSlackAndDual(const Robot& robot, const double dtau, 
                          const Eigen::MatrixXd& du_dq,
                          const Eigen::MatrixXd& du_dv,
                          const Eigen::MatrixXd& du_da,
                          const Eigen::VectorXd& dq,
                          const Eigen::VectorXd& dv,
                          const Eigen::VectorXd& du); 

  void condenseSlackAndDual(const Robot& robot, const double dtau, 
                            const Eigen::MatrixXd& du_dq,
                            const Eigen::MatrixXd& du_dv,
                            const Eigen::MatrixXd& du_da, Eigen::MatrixXd& Cqq, 
                            Eigen::MatrixXd& Cqv, Eigen::MatrixXd& Cqa, 
                            Eigen::MatrixXd& Cvv, Eigen::MatrixXd& Cva, 
                            Eigen::MatrixXd& Caa, Eigen::VectorXd& Cq,
                            Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void augmentDualResidual(const Robot& robot, const double dtau, 
                           const Eigen::MatrixXd& du_dq,
                           const Eigen::MatrixXd& du_dv, 
                           const Eigen::MatrixXd& du_da, Eigen::VectorXd& Cq,
                           Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  double squaredConstraintsResidualNrom(const Robot& robot, const double dtau, 
                                        const Eigen::VectorXd& u);

private:
  unsigned int dimq_, dimv_, dimc_;
  double barrier_;
  Eigen::VectorXd umin_, slack_, dual_, residual_, FB_residual_, dslack_, 
                  ddual_, newton_residual_;
  Eigen::MatrixXd partial_du_;
};

} // namespace idocp


#endif // IDOCP_JOINT_TORQUE_LOWER_LIMITS_HPP_