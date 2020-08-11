#ifndef IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_SPACE_CONSTRAINTS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/constraints/joint_space_constraints/joint_variables_upper_limits.hpp"
#include "idocp/constraints/joint_space_constraints/joint_variables_lower_limits.hpp"


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

  // Use default move constructor.
  JointSpaceConstraints(JointSpaceConstraints&&) noexcept = default;

  // Use default move assign operator.
  JointSpaceConstraints& operator=(JointSpaceConstraints&&) noexcept = default;

  void setTimeStep(const int time_step);

  bool isFeasible(const Eigen::Ref<const Eigen::VectorXd>& q, 
                  const Eigen::Ref<const Eigen::VectorXd>& v, 
                  const Eigen::Ref<const Eigen::VectorXd>& a, 
                  const Eigen::Ref<const Eigen::VectorXd>& u);

  void setSlackAndDual(const double dtau, 
                       const Eigen::Ref<const Eigen::VectorXd>& q, 
                       const Eigen::Ref<const Eigen::VectorXd>& v, 
                       const Eigen::Ref<const Eigen::VectorXd>& a, 
                       const Eigen::Ref<const Eigen::VectorXd>& u);

  void augmentDualResidual(const double dtau, Eigen::Ref<Eigen::VectorXd> Cu);

  void augmentDualResidual(const double dtau, Eigen::Ref<Eigen::VectorXd> Cq, 
                           Eigen::Ref<Eigen::VectorXd> Cv, 
                           Eigen::Ref<Eigen::VectorXd> Ca);

  void condenseSlackAndDual(const double dtau, 
                            const Eigen::Ref<const Eigen::VectorXd>& u, 
                            Eigen::Ref<Eigen::MatrixXd> Cuu, 
                            Eigen::Ref<Eigen::VectorXd> Cu);

  void condenseSlackAndDual(const double dtau, 
                            const Eigen::Ref<const Eigen::VectorXd>& q, 
                            const Eigen::Ref<const Eigen::VectorXd>& v, 
                            const Eigen::Ref<const Eigen::VectorXd>& a, 
                            Eigen::Ref<Eigen::MatrixXd> Cqq, 
                            Eigen::Ref<Eigen::MatrixXd> Cvv, 
                            Eigen::Ref<Eigen::MatrixXd> Caa,  
                            Eigen::Ref<Eigen::VectorXd> Cq, 
                            Eigen::Ref<Eigen::VectorXd> Cv, 
                            Eigen::Ref<Eigen::VectorXd> Ca);

  void computeSlackAndDualDirection(
      const double dtau, const Eigen::Ref<const Eigen::VectorXd>& dq, 
      const Eigen::Ref<const Eigen::VectorXd>& dv, 
      const Eigen::Ref<const Eigen::VectorXd>& da, 
      const Eigen::Ref<const Eigen::VectorXd>& du);

  double maxSlackStepSize();

  double maxDualStepSize();

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  double residualL1Nrom(const double dtau, 
                        const Eigen::Ref<const Eigen::VectorXd>& q, 
                        const Eigen::Ref<const Eigen::VectorXd>& v, 
                        const Eigen::Ref<const Eigen::VectorXd>& a, 
                        const Eigen::Ref<const Eigen::VectorXd>& u);

  double residualSquaredNrom(const double dtau, 
                             const Eigen::Ref<const Eigen::VectorXd>& q, 
                             const Eigen::Ref<const Eigen::VectorXd>& v, 
                             const Eigen::Ref<const Eigen::VectorXd>& a, 
                             const Eigen::Ref<const Eigen::VectorXd>& u);

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