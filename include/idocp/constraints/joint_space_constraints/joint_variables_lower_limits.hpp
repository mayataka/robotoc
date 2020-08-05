#ifndef IDOCP_CONSTRAINTS_JOINT_VARIABLES_LOWER_LIMITS_HPP_
#define IDOCP_CONSTRAINTS_JOINT_VARIABLES_LOWER_LIMITS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {
namespace pdipm {

class JointVariablesLowerLimits {
public:
  JointVariablesLowerLimits(const Robot& robot, const Eigen::VectorXd& xmin, 
                            const double barrier);

  JointVariablesLowerLimits();

  ~JointVariablesLowerLimits();

  // Use default copy constructor.
  JointVariablesLowerLimits(const JointVariablesLowerLimits&) = default;

  // Use default copy operator.
  JointVariablesLowerLimits& operator=(const JointVariablesLowerLimits&) 
      = default;

  // Use default move constructor.
  JointVariablesLowerLimits(JointVariablesLowerLimits&&) noexcept = default;

  // Use default move assign operator.
  JointVariablesLowerLimits& operator=(JointVariablesLowerLimits&&) noexcept 
      = default;

  bool isFeasible(const Eigen::VectorXd& v);

  void setSlackAndDual(const double dtau, const Eigen::VectorXd& x);

  void condenseSlackAndDual(const double dtau, const Eigen::VectorXd& x, 
                            Eigen::MatrixXd& Cxx, Eigen::VectorXd& Cx);

  void computeSlackAndDualDirection(const double dtau,
                                    const Eigen::VectorXd& dx);

  double maxSlackStepSize(const double margin_rate);

  double maxDualStepSize(const double margin_rate);

  void updateSlack(const double step_size);

  void updateDual(const double step_size);

  double costSlackBarrier();

  double costSlackBarrier(const double step_size);

  void augmentDualResidual(const double dtau, Eigen::VectorXd& Cx);

  double residualL1Nrom(const double dtau, const Eigen::VectorXd& x);

  double residualSquaredNrom(const double dtau, const Eigen::VectorXd& x);

private:
  int dimv_, dimc_, dim_passive_;
  bool has_floating_base_;
  double barrier_;
  Eigen::VectorXd xmin_, slack_, dual_, residual_, duality_, dslack_, ddual_;
};

} // namespace pdipm
} // namespace idocp


#endif // IDOCP_CONSTRAINTS_JOINT_VARIABLES_LOWER_LIMITS_HPP_