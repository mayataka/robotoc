#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class KKTMatrix
/// @brief The KKT matrix of a time stage.
///
class KKTMatrix {
public:
  ///
  /// @brief Construct a KKT matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  KKTMatrix(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  KKTMatrix();

  ///
  /// @brief Destructor. 
  ///
  ~KKTMatrix();

  ///
  /// @brief Default copy constructor. 
  ///
  KKTMatrix(const KKTMatrix&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  KKTMatrix& operator=(const KKTMatrix&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  KKTMatrix(KKTMatrix&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  KKTMatrix& operator=(KKTMatrix&&) noexcept = default;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  void setContactStatus(const Robot& robot);

  ///
  /// @brief Jacobian of the equality constraint with respect to acceleration a.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimc() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Ca();

  ///
  /// @brief Jacobian of the equality constraint with respect to the stack of 
  /// contact forces f.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimc() x KKTResidual::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Cf();

  ///
  /// @brief Jacobian of the equality constraint with respect to configuration q.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimc() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq();

  ///
  /// @brief Jacobian of the equality constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimc() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cv();

  ///
  /// @brief Jacobian of the floating base constraint with respect to 
  /// acceleration a.
  /// @return Reference to the Jacobian. Size is 
  /// Robot::dim_passive() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Ca_floating_base();

  ///
  /// @brief Jacobian of the floating base constraint with respect to the  
  /// stack of contact forces f.
  /// @return Reference to the Jacobian. Size is 
  /// Robot::dim_passive() x KKTResidual::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Cf_floating_base();

  ///
  /// @brief Jacobian of the floating base constraint with respect to 
  /// configuration q.
  /// @return Reference to the Jacobian. Size is 
  /// Robot::dim_passive() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq_floating_base();

  ///
  /// @brief Jacobian of the floating base constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// Robot::dim_passive() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cv_floating_base();

  ///
  /// @brief Jacobian of the contact constraint with respect to acceleration a.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Ca_contacts();

  ///
  /// @brief Jacobian of the contact constraint with respect to the stack of 
  /// contact forces f.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimf() x KKTResidual::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Cf_contacts();

  ///
  /// @brief Jacobian of the contact constraint with respect to configuration q.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq_contacts();

  ///
  /// @brief Jacobian of the contact constraint with respect to velocity v.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cv_contacts();

  ///
  /// @brief Jacobian of the equality constraint with respect to acceleration 
  /// and the stack of contact force, a and f.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimc() x (Robot::dimv()+KKTMatrix::dimf()).
  ///
  Eigen::Block<Eigen::MatrixXd> Caf();

  ///
  /// @brief Jacobian of the equality constraint with respect to configuration 
  /// and velocity q and v.
  /// @return Reference to the Jacobian. Size is 
  /// KKTMatrix::dimc() x 2*Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cqv();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration a. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaa();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and the 
  /// stack of the contact forces, a and f. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x KKTMatrix::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaf();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaq();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qav();

  ///
  /// @brief Hessian of the Lagrangian with respect to the stack of contact 
  /// forces and acceleration, f and a.
  /// @return Reference to the Hessian. Size is 
  /// KKTMatrix::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qfa();

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
  /// acceleration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqa();

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
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qva();

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
  /// 2*Robot::dimv() x 2*Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxx();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and the 
  /// stack of contact forces, (a, v) and (a, f). 
  /// @return Reference to the Hessian. Size is 
  /// (Robot::dimv()+KKTMatrix::dimf()) x (Robot::dimv()+KKTMatrix::dimf()).
  ///
  Eigen::Block<Eigen::MatrixXd> Qafaf();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and contact 
  /// forces, and configuraiton and velocity, (a, v) and (q, v). 
  /// @return Reference to the Hessian. Size is 
  /// (Robot::dimv()+KKTMatrix::dimf()) x 2*Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qafqv();

  ///
  /// @brief Hessian of the Lagrangian. 
  /// @return Reference to the Hessian. Size is 
  /// (3*Robot::dimv()+KKTMatrix::dimf()) x (3*Robot::dimv()+KKTMatrix::dimf()).
  ///
  Eigen::Block<Eigen::MatrixXd> costHessian();

  ///
  /// @brief Jacobian of the equality constraint. 
  /// @return Reference to the Hessian. Size is 
  /// KKTMatrix::dimc() x (3*Robot::dimv()+KKTMatrix::dimf()).
  ///
  Eigen::Block<Eigen::MatrixXd> constraintsJacobian();

  ///
  /// @brief Symmetrize the Hessian for matrix inversion. 
  ///
  void symmetrize();

  ///
  /// @brief Invert the KKT matrix. 
  /// @param[in] dtau Time step of the horzion.
  /// @param[out] kkt_matrix_inverse Inverse of the KKT matrix. Size must 
  /// be KKTMatrix::dimKKT() x KKTMatrix::dimKKT().
  ///
  template <typename MatrixType>
  void invert(const double dtau, 
              const Eigen::MatrixBase<MatrixType>& kkt_matrix_inverse);

  ///
  /// @brief Set the KKT residual zero.
  ///
  void setZeroMinimum();

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

  /// @brief Hessian of the Lagrangian with respect to the control input 
  //// torques u.
  Eigen::MatrixXd Quu;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration q.
  Eigen::MatrixXd Fqq;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration of the previous time step q_prev.
  Eigen::MatrixXd Fqq_prev;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd C_, Q_, Sc_, Sx_, FMinv_, C_H_inv_;
  bool has_floating_base_;
  int dimv_, dimx_, dim_passive_, dimf_, dimc_, a_begin_, f_begin_, q_begin_, 
      v_begin_, dimQ_, max_dimKKT_;
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

#include "idocp/ocp/kkt_matrix.hxx"

#endif // IDOCP_KKT_MATRIX_HPP_