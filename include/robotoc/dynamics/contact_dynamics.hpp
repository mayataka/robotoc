#ifndef ROBOTOC_CONTACT_DYNAMICS_HPP_
#define ROBOTOC_CONTACT_DYNAMICS_HPP_

#include <limits>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_solution.hpp"
#include "robotoc/core/split_direction.hpp"
#include "robotoc/core/split_kkt_residual.hpp"
#include "robotoc/core/split_kkt_matrix.hpp"
#include "robotoc/core/switching_constraint_residual.hpp"
#include "robotoc/core/switching_constraint_jacobian.hpp"
#include "robotoc/dynamics/contact_dynamics_data.hpp"


namespace robotoc {

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
  /// @brief Default destructor. 
  ///
  ~ContactDynamics() = default;

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
  static void eval(Robot& robot, const ContactStatus& contact_status,
                   ContactDynamicsData& data, const SplitSolution& s);

  ///
  /// @brief Computes the residual and derivatives of the contact dynamics  
  /// constraint and derivatives of it. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  static void linearize(Robot& robot, const ContactStatus& contact_status, 
                        ContactDynamicsData& data, const SplitSolution& s, 
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
  static void condense(Robot& robot, const ContactStatus& contact_status, 
                       ContactDynamicsData& data, const double dt, 
                       SplitKKTMatrix& kkt_matrix, 
                       SplitKKTResidual& kkt_residual);

  ///
  /// @brief Condenses the switching constraint. 
  /// @param[in, out] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] sc_residual Residual of the switching constraint. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  ///
  static void condense(const ContactDynamicsData& data, 
                       SwitchingConstraintJacobian& sc_jacobian,
                       SwitchingConstraintResidual& sc_residual,
                       SplitKKTMatrix& kkt_matrix);

  ///
  /// @brief Expands the primal variables, i.e., computes the Newton direction 
  /// of the condensed primal variables (acceleration a and the contact forces 
  /// f) of this stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  static void expandPrimal(const ContactDynamicsData& data, SplitDirection& d);

  ///
  /// @brief Expands the dual variables, i.e., computes the Newton direction 
  /// of the condensed dual variables (Lagrange multipliers) of this stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] dts Direction of the switching time regarding of this time stage. 
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  static void expandDual(ContactDynamicsData& data, 
                         const double dt, const double dts, 
                         const SplitDirection& d_next, SplitDirection& d);

  ///
  /// @brief Expands the dual variables, i.e., computes the Newton direction 
  /// of the condensed dual variables (Lagrange multipliers) of this stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] dts Direction of the switching time regarding of this time stage. 
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] d Split direction of this time stage.
  /// 
  static void expandDual(ContactDynamicsData& data, 
                         const double dt, const double dts, 
                         const SplitDirection& d_next, 
                         const SwitchingConstraintJacobian& sc_jacobian,
                         SplitDirection& d);

  ///
  /// @brief Computes the residual in the contact dynamics constraint. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] s Split solution of this time stage.
  ///
  void evalContactDynamics(Robot& robot, const ContactStatus& contact_status,
                           const SplitSolution& s);

  ///
  /// @brief Computes the residual and derivatives of the contact dynamics  
  /// constraint and derivatives of it. 
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status of this time stage. 
  /// @param[in] s Split solution of this time stage.
  /// @param[in, out] kkt_residual Split KKT residual of this time stage.
  ///
  void linearizeContactDynamics(Robot& robot, 
                                const ContactStatus& contact_status, 
                                const SplitSolution& s, 
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
  /// @brief Condenses the switching constraint. 
  /// @param[in, out] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] sc_residual Residual of the switching constraint. 
  /// @param[in, out] kkt_matrix Split KKT matrix of this time stage.
  ///
  void condenseSwitchingConstraint(SwitchingConstraintJacobian& sc_jacobian,
                                   SwitchingConstraintResidual& sc_residual,
                                   SplitKKTMatrix& kkt_matrix) const;

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
  /// @param[in] dts Direction of the switching time regarding of this time stage. 
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void expandDual(const double dt, const double dts, 
                  const SplitDirection& d_next, SplitDirection& d);

  ///
  /// @brief Expands the dual variables, i.e., computes the Newton direction 
  /// of the condensed dual variables (Lagrange multipliers) of this stage.
  /// @param[in] dt Time step of this time stage. 
  /// @param[in] dts Direction of the switching time regarding of this time stage. 
  /// @param[in] d_next Split direction of the next stage.
  /// @param[in] sc_jacobian Jacobian of the switching constraint. 
  /// @param[in, out] d Split direction of this time stage.
  /// 
  void expandDual(const double dt, const double dts, 
                  const SplitDirection& d_next, 
                  const SwitchingConstraintJacobian& sc_jacobian,
                  SplitDirection& d);

  ///
  /// @brief Returns the squared norm of the KKT residual, that is, 
  /// the primal and dual residual of the contact dynamics constraint. 
  /// @return Squared norm of the KKT residual in the contact dynamics 
  /// constraint.
  ///
  double KKTError() const {
    return (data_.IDC().squaredNorm() + data_.lu_passive.squaredNorm());
  }

  ///
  /// @brief Returns the lp norm of the constraint violation, that is,
  /// the primal residual in the contact dynamics. Default norm is l1-norm.
  /// You can specify l-infty norm by passing Eigen::Infinity as the 
  /// template parameter.
  /// @tparam p Index of norm. Default is 1 (l1-norm).
  /// @return The lp norm of the constraint violation.
  ///
  template <int p=1>
  double constraintViolation() const {
    return data_.IDC().template lpNorm<p>();
  }

private:
  ContactDynamicsData data_;
  static constexpr int kDimFloatingBase = 6;

};

} // namespace robotoc 

#endif // ROBOTOC_CONTACT_DYNAMICS_HPP_ 