#ifndef IDOCP_FORWARD_SWITCHING_CONSTRAINT_HPP_ 
#define IDOCP_FORWARD_SWITCHING_CONSTRAINT_HPP_

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_state_constraint_jacobian.hpp"


namespace idocp {

class ForwardSwitchingConstraint {
public:
  ForwardSwitchingConstraint(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ForwardSwitchingConstraint();

  ///
  /// @brief Destructor. 
  ///
  ~ForwardSwitchingConstraint();

  ///
  /// @brief Default copy constructor. 
  ///
  ForwardSwitchingConstraint(const ForwardSwitchingConstraint&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  ForwardSwitchingConstraint& operator=(
      const ForwardSwitchingConstraint&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ForwardSwitchingConstraint(ForwardSwitchingConstraint&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ForwardSwitchingConstraint& operator=(
      ForwardSwitchingConstraint&&) noexcept = default;

  void linearizeSwitchingConstraint(Robot& robot, 
                                    const ImpulseStatus& impulse_status, 
                                    const double dt1, const double dt2, 
                                    const SplitSolution& s, 
                                    SplitKKTMatrix& kkt_matrix, 
                                    SplitKKTResidual& kkt_residual, 
                                    SplitStateConstraintJacobian& jac);

  void computeSwitchingConstraintResidual(Robot& robot, 
                                          const ImpulseStatus& impulse_status, 
                                          const double dt1, const double dt2, 
                                          const SplitSolution& s, 
                                          SplitKKTResidual& kkt_residual);

  static double squaredNormSwitchingConstraintResidual(
      const SplitKKTResidual& kkt_residual);

  static double l1NormSwitchingConstraintResidual(
      const SplitKKTResidual& kkt_residual);

private:
  Eigen::VectorXd q_, dq_;

};

} // namespace idocp

#include "idocp/ocp/forward_switching_constraint.hxx"

#endif // IDOCP_FORWARD_SWITCHING_CONSTRAINT_HPP_ 