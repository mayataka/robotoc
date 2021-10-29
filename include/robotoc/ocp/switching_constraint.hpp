#ifndef ROBOTOC_SWITCHING_CONSTRAINT_HPP_ 
#define ROBOTOC_SWITCHING_CONSTRAINT_HPP_

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"
#include "robotoc/ocp/split_solution.hpp"
#include "robotoc/ocp/split_kkt_matrix.hpp"
#include "robotoc/ocp/split_kkt_residual.hpp"
#include "robotoc/ocp/switching_constraint_residual.hpp"
#include "robotoc/ocp/switching_constraint_jacobian.hpp"


namespace robotoc {

///
/// @class SwitchingConstraint
/// @brief The pure-state constraint representing switching constraint that is 
///  transformed into the mixed state-control equality constraint .
///
class SwitchingConstraint {
public:
  ///
  /// @brief Constructs a switching constraint.
  /// @param[in] robot Robot model. 
  ///
  SwitchingConstraint(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  SwitchingConstraint();

  ///
  /// @brief Destructor. 
  ///
  ~SwitchingConstraint();

  ///
  /// @brief Default copy constructor. 
  ///
  SwitchingConstraint(const SwitchingConstraint&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SwitchingConstraint& operator=(const SwitchingConstraint&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SwitchingConstraint(SwitchingConstraint&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SwitchingConstraint& operator=(SwitchingConstraint&&) noexcept = default;

  ///
  /// @brief Computes the residual in the switching constraint, i.e., the 
  /// contact position constraint. Note that the internal kinematics data 
  /// of robot is updated.  
  /// @param[in] robot Robot model. Kinematics must be updated.
  /// @param[in] impulse_status Impulse status. 
  /// @param[in] dt1 Time step of the time stage 2 stage before the impulse.
  /// @param[in] dt2 Time step of the time stage just before the impulse.
  /// @param[in] s Split solution of the time stage 2 stage before the impulse.
  /// @param[in, out] sc_residual Residual of the switching constraint. 
  ///
  void evalSwitchingConstraint(Robot& robot, const ImpulseStatus& impulse_status, 
                               const double dt1, const double dt2, 
                               const SplitSolution& s, 
                               SwitchingConstraintResidual& sc_residual);

  ///
  /// @brief Linearizes the switching constraint, i.e., the contact position 
  /// constraint. Note that the internal kinematics data of robot is updated.
  /// @param[in] robot Robot model. Kinematics must be updated.
  /// @param[in] impulse_status Impulse status. 
  /// @param[in] dt1 Time step of the time stage 2 stage before the impulse.
  /// @param[in] dt2 Time step of the time stage just before the impulse.
  /// @param[in] s Split solution of the time stage 2 stage before the impulse.
  /// @param[in, out] kkt_matrix Split KKT matrix of the time stage 2 stage
  /// before the impulse.
  /// @param[in, out] kkt_residual Split KKT residual of the time stage 2 stage
  /// before the impulse.
  /// @param[in, out] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] sc_residual Residual of the switching constraint. 
  ///
  void linearizeSwitchingConstraint(
      Robot& robot, const ImpulseStatus& impulse_status, const double dt1, 
      const double dt2, const SplitSolution& s, SplitKKTMatrix& kkt_matrix, 
      SplitKKTResidual& kkt_residual, 
      SwitchingConstraintJacobian& sc_jacobian,
      SwitchingConstraintResidual& sc_residual);

private:
  Eigen::VectorXd q_, dq_;
  bool has_floating_base_;

};

} // namespace robotoc

#include "robotoc/ocp/switching_constraint.hxx"

#endif // ROBOTOC_SWITCHING_CONSTRAINT_HPP_ 