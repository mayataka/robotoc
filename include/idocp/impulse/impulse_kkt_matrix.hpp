#ifndef IDOCP_IMPULSE_KKT_MATRIX_HPP_
#define IDOCP_IMPULSE_KKT_MATRIX_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/schur_complement.hpp"


namespace idocp {

///
/// @class ImpulseKKTMatrix
/// @brief The KKT matrix of a impulse time stage.
///
class ImpulseKKTMatrix {
public:
  ///
  /// @brief Construct a impulse KKT matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] use_contact_position_constraint true if you treat the contact 
  /// position constraint in this impulse stage. false if not.
  ///
  ImpulseKKTMatrix(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseKKTMatrix();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseKKTMatrix();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseKKTMatrix(const ImpulseKKTMatrix&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseKKTMatrix& operator=(const ImpulseKKTMatrix&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  ImpulseKKTMatrix(ImpulseKKTMatrix&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseKKTMatrix& operator=(ImpulseKKTMatrix&&) noexcept = default;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Jacobian of the contact position and velocity constraint with  
  /// respect to configuration q.
  /// @return Reference to the Jacobian. Size is 
  /// 2 * ContactStatus::dimf() x KKTMatrix::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq();

  ///
  /// @brief Jacobian of the contact position and velocity constraint with  
  /// respect to generalized velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// 2 * ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cv();

  ///
  /// @brief Jacobian of the contact position and velocity constraint with  
  /// respect to q and v.
  /// @return Reference to the Jacobian. Size is 
  /// 2 * ContactStatus::dimf() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cqv();

  ///
  /// @brief Jacobian of the contact position and velocity constraint with  
  /// respect to generalized velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// 2 * ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cf();

  ///
  /// @brief Jacobian of the contact position constraint with respect to 
  /// configuration q.
  /// @return Reference to the Jacobian. Size is 
  /// ContactStatus::dimf() x KKTMatrix::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq_contact_position();

  ///
  /// @brief Jacobian of the contact position constraint with respect to 
  /// generalized velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cv_contact_position();

  ///
  /// @brief Jacobian of the contact position constraint with respect to 
  /// generalized velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cf_contact_position();

  ///
  /// @brief Jacobian of the contact velocity constraint with respect to 
  /// configuration q.
  /// @return Reference to the Jacobian. Size is 
  /// ContactStatus::dimf() x KKTMatrix::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq_contact_velocity();

  ///
  /// @brief Jacobian of the contact velocity constraint with respect to 
  /// generalized velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cv_contact_velocity();

  ///
  /// @brief Jacobian of the contact velocity constraint with respect to 
  /// generalized velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cf_contact_velocity();

  ///
  /// @brief Hessian of the Lagrangian with respect to the stack of contact 
  /// forces f.
  /// @return Reference to the Hessian. Size is 
  /// KKTMatrix::dimf() x KKTMatrix::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qff();

  ///
  /// @brief Hessian of the Lagrangian with respect to the stack of contact 
  /// forces and configuration, f and q.
  /// @return Reference to the Hessian. Size is 
  /// KKTMatrix::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qfq();

  ///
  /// @brief Hessian of the Lagrangian with respect to the stack of contact 
  /// forces and velocity, f and v.
  /// @return Reference to the Hessian. Size is 
  /// KKTMatrix::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qfv();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and 
  /// the stack of contact forces, q and f. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x KKTMatrix::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqf();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqq();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and 
  /// the stack of contact forces, v and f. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x KKTMatrix::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvf();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and 
  /// configuration, v and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxx();

  ///
  /// @brief Hessian of the Lagrangian. 
  /// @return Reference to the Hessian. Size is 
  /// (3 * Robot::dimv() + KKTMatrix::dimf()) x (3 * Robot::dimv() + KKTMatrix::dimf()).
  ///
  Eigen::Block<Eigen::MatrixXd> costHessian();

  ///
  /// @brief Jacobian of the equality constraint. 
  /// @return Reference to the Hessian. Size is 
  /// KKTMatrix::dimc() x (3 * Robot::dimv() + KKTMatrix::dimf()).
  ///
  Eigen::Block<Eigen::MatrixXd> constraintsJacobian();

  ///
  /// @brief Symmetrize the Hessian for matrix inversion. 
  ///
  void symmetrize();

  ///
  /// @brief Invert the condensed KKT matrix. 
  /// @param[out] kkt_matrix_inverse Inverse of the KKT matrix. Size must 
  /// be KKTMatrix::dimKKT() x KKTMatrix::dimKKT().
  ///
  template <typename MatrixType>
  void invert(const Eigen::MatrixBase<MatrixType>& kkt_matrix_inverse);

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

  /// @brief Hessian of the Lagrangian with respect to the impulse chane in 
  /// the generalized velocity dv.
  Eigen::MatrixXd Qdvdv;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration q.
  Eigen::MatrixXd Fqq;

  /// @brief Derivative of the state equation with respect to the generalized 
  /// velocity v and configuration q.
  Eigen::MatrixXd Fvq;

  /// @brief Derivative of the state equation with respect to the generalized 
  /// velocity v.
  Eigen::MatrixXd Fvv;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration of the previous time step q_prev.
  Eigen::MatrixXd Fqq_prev;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SchurComplement schur_complement_;
  Eigen::MatrixXd C_, Q_, Sx_, FMinv_;
  bool has_floating_base_;
  int dimv_, dimx_, dimf_, dimc_, f_begin_, q_begin_, v_begin_, dimQ_, 
      dimKKT_, max_dimKKT_;
  static constexpr int kDimFloatingBase = 6;

  ///
  /// @brief Invert the cost Hessian matrix. 
  /// @param[out] hessian_inverse Inverse of the Hessian matrix.
  ///
  template <typename MatrixType>
  void invertConstrainedHessian(
      const Eigen::MatrixBase<MatrixType>& hessian_inverse);

};

} // namespace idocp 

#include "idocp/impulse/impulse_kkt_matrix.hxx"

#endif // IDOCP_IMPULSE_KKT_MATRIX_HPP_ 