#ifndef IDOCP_JOINT_VELOCITY_LOWER_LIMITS_HPP_
#define IDOCP_JOINT_VELOCITY_LOWER_LIMITS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"


namespace idocp {

class JointVelocityLowerLimits {
public:
  JointVelocityLowerLimits(const Robot& robot);

  // Use default copy constructor.
  JointVelocityLowerLimits(const JointVelocityLowerLimits&) = default;

  // Use default copy operator.
  JointVelocityLowerLimits& operator=(const JointVelocityLowerLimits&) = default;

  void setSlackAndDual(const Robot& robot, const double barrier, 
                       const Eigen::VectorXd& v);

  void linearizeConstraint(const Robot& robot, const double barrier, 
                           const Eigen::VectorXd& v);

  void updateSlackAndDual(const Robot& robot, const Eigen::VectorXd& dv); 

  void condenseSlackAndDual(const Robot& robot, Eigen::MatrixXd& Cvv,     
                            Eigen::VectorXd& Cv) const;

  void augmentDualResidual(const Robot& robot, Eigen::VectorXd& Cv);

  double squaredConstraintsResidualNrom(const Robot& robot, 
                                        const Eigen::VectorXd& v);

private:
  unsigned int dimq_, dimv_, dimc_;
  Eigen::VectorXd vmin_, slack_, dual_, residual_, FB_residual_, dslack_, 
                  ddual_;
};

} // namespace idocp


#endif // IDOCP_JOINT_VELOCITY_LOWER_LIMITS_HPP_