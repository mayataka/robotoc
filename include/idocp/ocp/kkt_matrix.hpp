#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/ocp/schur_complement.hpp"


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
  /// @brief Set contact status, i.e., set dimension of the contacts.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] impulse_status Contact status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Jacobian of the state equation of the configuration with respect  
  /// to the control input torques. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqu();

  ///
  /// @brief const version of KKTMatrix::Fqu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqu() const;

  ///
  /// @brief Jacobian of the state equation of the configuration with respect  
  /// to the configuration. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqq();

  ///
  /// @brief const version of KKTMatrix::Fqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqq() const;

  ///
  /// @brief Jacobian of the state equation of the configuration with respect  
  /// to the generalized velocity. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqv();

  ///
  /// @brief const version of KKTMatrix::Fqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqv() const;

  ///
  /// @brief Jacobian of the state equation of the generalized velocity with 
  /// respect to the control input torques. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvu();

  ///
  /// @brief const version of KKTMatrix::Fvu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvu() const;

  ///
  /// @brief Jacobian of the state equation of the generalized velocity with 
  /// respect to the configuration. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvq();

  ///
  /// @brief const version of KKTMatrix::Fvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvq() const;

  ///
  /// @brief Jacobian of the state equation of the generalized velocity with 
  /// respect to the generalized velocity. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvv();

  ///
  /// @brief const version of KKTMatrix::Fvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvv() const;

  ///
  /// @brief Jacobian of the state equation with respect to the control input 
  /// torques.
  /// @return Reference to the block part of the Hessian. 
  /// Size is 2 * Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Fxu();

  ///
  /// @brief const version of KKTMatrix::Fxu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fxu() const;

  ///
  /// @brief Jacobian of the state equation with respect to the state, i.e.,
  /// the configuration and the generalized velocity.
  /// @return Reference to the block part of the Hessian. 
  /// Size is 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fxx();

  ///
  /// @brief const version of KKTMatrix::Fxx().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fxx() const;

  ///
  /// @brief Jacobian of the contact position constraint related to impulse 
  /// condition with respect to the configuration. KKTMatrix::setImpulseStatus()
  /// muset be called to set the impulse dimension before calling this function.
  /// @return Reference to the block part of the Hessian. 
  /// Size is ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Cq();

  ///
  /// @brief const version of KKTMatrix::Cq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Cq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_full();

  ///
  /// @brief const version of KKTMatrix::Quu_full().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dim_passive() x Robot::dim_passive().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_passive_topLeft();

  ///
  /// @brief const version of KKTMatrix::Quu_passive_topLeft().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_passive_topLeft() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimu() x Robot::dim_passive().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_passive_topRight();

  ///
  /// @brief const version of KKTMatrix::Quu_passive_topRight().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_passive_topRight() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dim_passive() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_passive_bottomLeft();

  ///
  /// @brief const version of KKTMatrix::Quu_passive_bottomLeft().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_passive_bottomLeft() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu();

  ///
  /// @brief const version of KKTMatrix::Quu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq_full();

  ///
  /// @brief const version of KKTMatrix::Quq_full().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dim_passive() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq_passive();

  ///
  /// @brief const version of KKTMatrix::Quq_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq();

  ///
  /// @brief const version of KKTMatrix::Quq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// generalized velocity. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv_full();

  ///
  /// @brief const version of KKTMatrix::Quv_full().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// generalized velocity. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dim_passive() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv_passive();

  ///
  /// @brief const version of KKTMatrix::Quv_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// generalized velocity. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv();

  ///
  /// @brief const version of KKTMatrix::Quv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and the 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu_full();

  ///
  /// @brief const version of KKTMatrix::Qqu_full().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and the 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x Robot::dim_passive().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu_passive();

  ///
  /// @brief const version of KKTMatrix::Qqu_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and the 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu();

  ///
  /// @brief const version of KKTMatrix::Qqu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqq();

  ///
  /// @brief const version of KKTMatrix::Qqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief const version of KKTMatrix::Qqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu_full();

  ///
  /// @brief const version of KKTMatrix::Qvu_full().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x Robot::dim_passive().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu_passive();

  ///
  /// @brief const version of KKTMatrix::Qvu_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu();

  ///
  /// @brief const version of KKTMatrix::Qvu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief const version of KKTMatrix::Qvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief const version of KKTMatrix::Qvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to control input torques
  /// and state. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qux_full();

  ///
  /// @brief const version of KKTMatrix::Qux().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qux_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to control input torques
  /// and state. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dim_passive() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qux_passive();

  ///
  /// @brief const version of KKTMatrix::Qux_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qux_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to control input torques
  /// and state. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimu() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qux();

  ///
  /// @brief const version of KKTMatrix::Qux().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qux() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state and control input
  /// torques. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxu_full();

  ///
  /// @brief const version of KKTMatrix::Qux_full().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state and control input
  /// torques. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x Robot::dim_passive().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxu_passive();

  ///
  /// @brief const version of KKTMatrix::Qxu_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state and control input
  /// torques. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxu();

  ///
  /// @brief const version of KKTMatrix::Qxu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxx();

  ///
  /// @brief const version of KKTMatrix::Qxx().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxx() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// contact forces. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() + ContactStatus::dimf() 
  /// x  Robot::dimv() + ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaaff();

  ///
  /// @brief const version of KKTMatrix::Qaaff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaaff() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaa();

  ///
  /// @brief const version of KKTMatrix::Qaa().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaa() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration. 
  /// @return Reference to the Hessian. Size is 
  /// ContactStatus::dimf() x ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qff();

  ///
  /// @brief const version of KKTMatrix::Qff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qff() const;

  ///
  /// @brief Symmetrize the Hessian for matrix inversion. 
  ///
  void symmetrize();

  ///
  /// @brief Invert the KKT matrix. 
  /// @param[out] KKT_matrix_inverse Inverse of the KKT matrix. Size must 
  /// be KKTMatrix::dimKKT() x KKTMatrix::dimKKT().
  ///
  template <typename MatrixType>
  void invert(const Eigen::MatrixBase<MatrixType>& KKT_matrix_inverse);

  ///
  /// @brief Set the all components zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the condensed KKT condition.
  /// @return Dimension of the condensed KKT condition.
  ///
  int dimKKT() const;

  ///
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of the stack of contact forces.
  ///
  int dimf() const;

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of the stack of impulse forces.
  ///
  int dimp() const;

  bool isApprox(const KKTMatrix& other) const;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration of the previous time step q_prev.
  Eigen::MatrixXd Fqq_prev;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SchurComplement schur_complement_;
  Eigen::MatrixXd F_, C_, Q_, Qaaff_full_;
  bool has_floating_base_;
  int dimv_, dimx_, dimu_, dim_passive_, dimf_, dimp_, u_begin_, q_begin_, 
      v_begin_, dimKKT_;
  static constexpr int kDimFloatingBase = 6;

};

} // namespace idocp 

#include "idocp/ocp/kkt_matrix.hxx"

#endif // IDOCP_KKT_MATRIX_HPP_