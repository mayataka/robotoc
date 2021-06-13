#ifndef IDOCP_CONTACT_DYNAMICS_HPP_
#define IDOCP_CONTACT_DYNAMICS_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_solution.hpp"
#include "idocp/ocp/split_direction.hpp"
#include "idocp/ocp/split_kkt_residual.hpp"
#include "idocp/ocp/split_kkt_matrix.hpp"
#include "idocp/ocp/contact_dynamics_data.hpp"
#include "idocp/ocp/split_switching_constraint_residual.hpp"
#include "idocp/ocp/split_switching_constraint_jacobian.hpp"


namespace idocp {

///
/// @class ContactDynamics
/// @brief Contact dynamics constraint.
///
class ContactDynamics {
public:
  ///
  /// @brief Constructs the contact dynamics.
  /// @param[in] robot Robot model. 
  ///
  ContactDynamics(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ContactDynamics();

  ///
  /// @brief Destructor. 
  ///
  ~ContactDynamics();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactDynamics(const ContactDynamics&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ContactDynamics& operator=(const ContactDynamics&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ContactDynamics(ContactDynamics&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactDynamics& operator=(ContactDynamics&&) noexcept = default;

  ///
  /// @brief Computes the residual in the contact dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] s Split solution of this time stage.
  ///
  void computeContactDynamicsResidual(Robot& robot, 
                                      const ContactStatus& contact_status,
                                      const SplitSolution& s);

  ///
  /// @brief Computes the residual and derivatives of the contact dynamics  
  /// constraint and derivatives of it. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void linearizeContactDynamics(Robot& robot, 
                                const ContactStatus& contact_status, 
                                const double dt, const SplitSolution& s, 
                                SplitKKTResidual& kkt_residual);

  ///
  /// @brief Condenses the acceleration, contact forces, and Lagrange
  /// multipliers. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] dt Time step of this time stage. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void condenseContactDynamics(Robot& robot, 
                               const ContactStatus& contact_status, 
                               const double dt, SplitKKTMatrix& kkt_matrix, 
                               SplitKKTResidual& kkt_residual);

  ///
  /// @brief Expands the primal variables, i.e., computes the Newton direction 
  /// of the condensed primal variables (acceleration a and the contact forces 
  /// f) of this stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void expandPrimal(SplitDirection& d) const;

  ///
  /// @brief Expands the dual variables, i.e., computes the Newton direction 
  /// of the condensed dual variables (Lagrange multipliers) of this stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  template <typename SplitDirectionType>
  void expandDual(const double dt, const SplitDirectionType& d_next, 
                  SplitDirection& d);

  ///
  /// @brief Condenses the switching constraint. 
  /// @param[in, out] switch_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] switch_residual Residual of the switching constraint. 
  ///
  void condenseSwitchingConstraint(
      SplitSwitchingConstraintJacobian& switch_jacobian,
      SplitSwitchingConstraintResidual& switch_residual) const;

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual of the contact dynamics constraint. 
  /// @param[in] dt Time step of this time stage. 
  /// @return Squared norm of the KKT residual in the contact dynamics 
  /// constraint.
  ///
  double squaredNormKKTResidual(const double dt) const;

  ///
  /// @brief Returns l1-norm of the constraint violation, that is, the primal
  /// residual in the contact dynamics constraint. 
  /// @return l1-norm of the constraint violation.
  ///
  double l1NormConstraintViolation() const;

private:
  ContactDynamicsData data_;
  bool has_floating_base_, has_active_contacts_;
  int dimv_, dimu_, dim_passive_;
  static constexpr int kDimFloatingBase = 6;

};

} // namespace idocp 

#include "idocp/ocp/contact_dynamics.hxx"

#endif // IDOCP_CONTACT_DYNAMICS_HPP_ 