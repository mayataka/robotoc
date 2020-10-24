#ifndef IDOCP_CONTACT_DYNAMICS_JACOBIAN_HPP_
#define IDOCP_CONTACT_DYNAMICS_JACOBIAN_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


namespace idocp {

///
/// @class ContactDynamicsJacobian
/// @brief The KKT matrix of a time stage.
///
class ContactDynamicsJacobian {
public:
  ///
  /// @brief Construct a KKT matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ContactDynamicsJacobian(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ContactDynamicsJacobian();

  ///
  /// @brief Destructor. 
  ///
  ~ContactDynamicsJacobian();

  ///
  /// @brief Default copy constructor. 
  ///
  ContactDynamicsJacobian(const ContactDynamicsJacobian&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ContactDynamicsJacobian& operator=(const ContactDynamicsJacobian&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  ContactDynamicsJacobian(ContactDynamicsJacobian&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ContactDynamicsJacobian& operator=(ContactDynamicsJacobian&&) noexcept = default;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Jacobian of the contact constraint with respect to acceleration a.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Ca_contacts();

  ///
  /// @brief Jacobian of the contact constraint with respect to the stack of 
  /// contact forces f.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x ContactDynamicsJacobian::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Cf_contacts();

  ///
  /// @brief Jacobian of the contact constraint with respect to configuration q.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq_contacts();

  ///
  /// @brief Jacobian of the contact constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cv_contacts();

  ///
  /// @brief Jacobian of the contact constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> dID_dq();

  ///
  /// @brief Jacobian of the contact constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> dID_dv();

  ///
  /// @brief Jacobian of the contact constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> dID_da();

  ///
  /// @brief Jacobian of the contact constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactDynamicsJacobian::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> dID_df();

  ///
  /// @brief Set the KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the KKT at the current contact status.
  /// @return Dimension of the KKT at the current contact status.
  ///
  int dimKKT() const;

  ///
  /// @brief Returns the maximum dimension of the KKT.
  /// @return Maximum dimension of the KKT at the current contact status.
  ///
  int max_dimKKT() const;

  ///
  /// @brief Returns the dimension of equality constraint at the current 
  /// contact status.
  /// @return Dimension of equality constraint.
  ///
  int dimc() const;

  ///
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of the stack of contact forces.
  ///
  int dimf() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SchurComplement schur_complement_;
  Eigen::MatrixXd C_, Q_, Sx_, FMinv_;
  bool has_floating_base_;
  int dimv_, dimx_, dim_passive_, dimf_, dimc_, a_begin_, f_begin_, q_begin_, 
      v_begin_, dimQ_, max_dimKKT_;

};

} // namespace idocp 

#include "idocp/ocp/contact_dynamics_jacobian.hxx"

#endif // IDOCP_CONTACT_DYNAMICS_JACOBIAN_HPP_ 