#ifndef IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "constraints/joint_space_constraints/joint_variables_upper_limits.hpp"
#include "constraints/joint_space_constraints/joint_variables_lower_limits.hpp"


namespace idocp {
namespace pdipm {

class JointSpaceConstraints {
public:
  JointSpaceConstraints(const Robot& robot);

  JointSpaceConstraints();

  ~JointSpaceConstraints();

  // Use default copy constructor.
  JointSpaceConstraints(const JointSpaceConstraints&) = default;

  // Use default copy operator.
  JointSpaceConstraints& operator=(const JointSpaceConstraints&) = default;

  void setTimeStep(const int time_step);

  bool isFeasible(const Eigen::VectorXd& q, const Eigen::VectorXd& v, 
                  const Eigen::VectorXd& a, const Eigen::VectorXd& u);

  void setSlackAndDual(const double dtau, const Eigen::VectorXd& q, 
                       const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                       const Eigen::VectorXd& u);

  void augmentDualResidual(const double dtau, Eigen::VectorXd& Cu);

  void augmentDualResidual(const double dtau, Eigen::VectorXd& Cq, 
                           Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void condenseSlackAndDual(const double dtau, const Eigen::VectorXd& q, 
                            const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                            Eigen::MatrixXd& Cqq, Eigen::MatrixXd& Cvv, 
                            Eigen::MatrixXd& Caa,  Eigen::VectorXd& Cq, 
                            Eigen::VectorXd& Cv, Eigen::VectorXd& Ca);

  void condenseSlackAndDual(const double dtau, const Eigen::VectorXd& u, 
                            Eigen::MatrixXd& Cuu, Eigen::VectorXd& Cu);

  void computeSlackAndDualDirection(const double dtau,
                                    const Eigen::VectorXd& dq,
                                    const Eigen::VectorXd& dv,
                                    const Eigen::VectorXd& da,
                                    const Eigen::VectorXd& du);

  double maxSlackStepSize();

  double maxDualStepSize();

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  double residualL1Nrom(const double dtau, const Eigen::VectorXd& q, 
                        const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                        const Eigen::VectorXd& u);

  double residualSquaredNrom(const double dtau, const Eigen::VectorXd& q, 
                             const Eigen::VectorXd& v, const Eigen::VectorXd& a, 
                             const Eigen::VectorXd& u);

private:
  int time_step_, dimq_, dimv_;
  double barrier_, fraction_to_boundary_margin_;
  pdipm::JointVariablesUpperLimits position_upper_limits_, 
                                   velocity_upper_limits_, 
                                   torque_upper_limits_;
  pdipm::JointVariablesLowerLimits position_lower_limits_, 
                                   velocity_lower_limits_, 
                                   torque_lower_limits_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_