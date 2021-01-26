#ifndef IDOCP_SPLIT_KKT_MATRIX_HPP_
#define IDOCP_SPLIT_KKT_MATRIX_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class SplitKKTMatrix
/// @brief The KKT matrix split into a time stage.
///
class SplitKKTMatrix {
public:
  using Matrix6d = Eigen::Matrix<double, 6, 6>;

  ///
  /// @brief Construct a KKT matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  SplitKKTMatrix(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitKKTMatrix();

  ///
  /// @brief Destructor. 
  ///
  ~SplitKKTMatrix();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitKKTMatrix(const SplitKKTMatrix&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitKKTMatrix& operator=(const SplitKKTMatrix&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  SplitKKTMatrix(SplitKKTMatrix&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitKKTMatrix& operator=(SplitKKTMatrix&&) noexcept = default;

  ///
  /// @brief Set contact status, i.e., set dimension of the contacts.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse, to zero.
  ///
  void setImpulseStatus();

  ///
  /// @brief Jacobian of the state equation of the configuration with respect  
  /// to the control input torques. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqu();

  ///
  /// @brief const version of SplitKKTMatrix::Fqu().
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
  /// @brief const version of SplitKKTMatrix::Fqq().
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
  /// @brief const version of SplitKKTMatrix::Fqv().
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
  /// @brief const version of SplitKKTMatrix::Fvu().
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
  /// @brief const version of SplitKKTMatrix::Fvq().
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
  /// @brief const version of SplitKKTMatrix::Fvv().
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
  /// @brief const version of SplitKKTMatrix::Fxu().
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
  /// @brief const version of SplitKKTMatrix::Fxx().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fxx() const;

  ///
  /// @brief Jacobian of the contact position constraint related to impulse 
  /// condition with respect to the configuration. 
  /// ImpulseSplitKKTMatrix::setImpulseStatus() muset be called to set the impulse 
  /// dimension before calling this function.
  /// @return Reference to the block part of the Hessian. 
  /// Size is ImpulseStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Pq();

  ///
  /// @brief const version of SplitKKTMatrix::Pq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Pq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_full();

  ///
  /// @brief const version of SplitKKTMatrix::Quu_full().
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
  /// @brief const version of SplitKKTMatrix::Quu_passive_topLeft().
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
  /// @brief const version of SplitKKTMatrix::Quu_passive_topRight().
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
  /// @brief const version of SplitKKTMatrix::Quu_passive_bottomLeft().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_passive_bottomLeft() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu();

  ///
  /// @brief const version of SplitKKTMatrix::Quu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq_full();

  ///
  /// @brief const version of SplitKKTMatrix::Quq_full().
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
  /// @brief const version of SplitKKTMatrix::Quq_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq();

  ///
  /// @brief const version of SplitKKTMatrix::Quq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// generalized velocity. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv_full();

  ///
  /// @brief const version of SplitKKTMatrix::Quv_full().
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
  /// @brief const version of SplitKKTMatrix::Quv_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input and 
  /// generalized velocity. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv();

  ///
  /// @brief const version of SplitKKTMatrix::Quv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and the 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu_full();

  ///
  /// @brief const version of SplitKKTMatrix::Qqu_full().
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
  /// @brief const version of SplitKKTMatrix::Qqu_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and the 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu();

  ///
  /// @brief const version of SplitKKTMatrix::Qqu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqq();

  ///
  /// @brief const version of SplitKKTMatrix::Qqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief const version of SplitKKTMatrix::Qqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu_full();

  ///
  /// @brief const version of SplitKKTMatrix::Qvu_full().
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
  /// @brief const version of SplitKKTMatrix::Qvu_passive().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// control input torques. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu();

  ///
  /// @brief const version of SplitKKTMatrix::Qvu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief const version of SplitKKTMatrix::Qvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief const version of SplitKKTMatrix::Qvv().
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
  /// @brief const version of SplitKKTMatrix::Qux().
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
  /// @brief const version of SplitKKTMatrix::Qux_passive().
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
  /// @brief const version of SplitKKTMatrix::Qux().
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
  /// @brief const version of SplitKKTMatrix::Qux_full().
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
  /// @brief const version of SplitKKTMatrix::Qxu_passive().
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
  /// @brief const version of SplitKKTMatrix::Qxu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxx();

  ///
  /// @brief const version of SplitKKTMatrix::Qxx().
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
  /// @brief const version of SplitKKTMatrix::Qaaff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaaff() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaa();

  ///
  /// @brief const version of SplitKKTMatrix::Qaa().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaa() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration. 
  /// @return Reference to the Hessian. Size is 
  /// ContactStatus::dimf() x ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qff();

  ///
  /// @brief const version of SplitKKTMatrix::Qff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qff() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to primal variables. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() + Robot::dimu() x 2 * Robot::dimv() + Robot::dimu().
  ///
  Eigen::Block<Eigen::MatrixXd> Qss();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Qss().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qss() const;

  ///
  /// @brief Constraint Jacobian. 
  /// @return Reference to the Jacobian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv() + Robot::dimu().
  ///
  Eigen::MatrixXd& Jac();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Jac().
  ///
  const Eigen::MatrixXd& Jac() const;


  ///
  /// @brief Symmetrize the Hessian for matrix inversion. 
  ///
  void symmetrize();

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
  int dimi() const;

  ///
  /// @brief Chech the equivalence of two SplitKKTMatrix.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitKKTMatrix& other) const;

  ///
  /// @brief Chech this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration of the previous time step q_prev.
  Eigen::MatrixXd Fqq_prev;

  /// @brief Inverse of the derivative of the state equation with respect to 
  /// the configuration. 
  Matrix6d Fqq_inv;

  /// @brief Inverse of the derivative of the state equation with respect to 
  /// the configuration of the previous time step q_prev.
  Matrix6d Fqq_prev_inv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::MatrixXd F_, Pq_full_, Q_, Qaaff_full_;
  bool has_floating_base_;
  int dimv_, dimx_, dimu_, dim_passive_, dimf_, dimi_, dimKKT_,
      u_begin_, q_begin_, v_begin_;
  static constexpr int kDimFloatingBase = 6;

};

} // namespace idocp 

#include "idocp/ocp/split_kkt_matrix.hxx"

#endif // IDOCP_SPLIT_KKT_MATRIX_HPP_ 