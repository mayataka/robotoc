#ifndef IDOCP_CONSTRAINTS_JOINT_VARIABLES_UPPER_LIMITS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_VARIABLES_UPPER_LIMITS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {
namespace pdipm {

class JointVariablesUpperLimits {
public:
  JointVariablesUpperLimits(const Robot& robot, 
                            const Eigen::Ref<const Eigen::VectorXd>& xmax, 
                            const double barrier);

  JointVariablesUpperLimits();

  ~JointVariablesUpperLimits();

  // Use default copy constructor.
  JointVariablesUpperLimits(const JointVariablesUpperLimits&) = default;

  // Use default copy operator.
  JointVariablesUpperLimits& operator=(const JointVariablesUpperLimits&) 
      = default;

  // Use default move constructor.
  JointVariablesUpperLimits(JointVariablesUpperLimits&&) noexcept = default;

  // Use default move assign operator.
  JointVariablesUpperLimits& operator=(JointVariablesUpperLimits&&) noexcept
      = default;

  bool isFeasible(const Eigen::Ref<const Eigen::VectorXd>& v);

  void setSlackAndDual(const double dtau, 
                       const Eigen::Ref<const Eigen::VectorXd>& x);

  void condenseSlackAndDual(const double dtau, 
                            const Eigen::Ref<const Eigen::VectorXd>& x, 
                            Eigen::Ref<Eigen::MatrixXd> Cxx, 
                            Eigen::Ref<Eigen::VectorXd> Cx);

  void computeSlackAndDualDirection(
      const double dtau, const Eigen::Ref<const Eigen::VectorXd>& dx);

  double maxSlackStepSize(const double margin_rate);

  double maxDualStepSize(const double margin_rate);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  void augmentDualResidual(const double dtau, Eigen::Ref<Eigen::VectorXd> Cx);

  double residualL1Nrom(const double dtau, 
                        const Eigen::Ref<const Eigen::VectorXd>& x);

  double residualSquaredNrom(const double dtau, 
                             const Eigen::Ref<const Eigen::VectorXd>& x);

private:
  int dimv_, dimc_, dim_passive_;
  bool has_floating_base_;
  double barrier_;
  Eigen::VectorXd xmax_, slack_, dual_, residual_, duality_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_JOINT_VARIABLES_UPPER_LIMITS_HPP_