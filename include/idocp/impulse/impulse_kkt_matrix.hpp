#ifndef IDOCP_IMPULSE_KKT_MATRIX_HPP_
#define IDOCP_IMPULSE_KKT_MATRIX_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"
#include "idocp/impulse/dynamic_schur_complement.hpp"


namespace idocp {

///
/// @class ImpulseKKTMatrix
/// @brief The KKT matrix of a time stage.
///
class ImpulseKKTMatrix {
public:
  ///
  /// @brief Construct a KKT matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
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
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] impulse_status Contact status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Jacobian of the state equation of the configuration with respect  
  /// to the impulse forces. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x ImpulseStatus::dimp().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqf();

  ///
  /// @brief const version of ImpulseKKTMatrix::Fqf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqf() const;

  ///
  /// @brief Jacobian of the state equation of the configuration with respect  
  /// to the configuration. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqq();

  ///
  /// @brief const version of ImpulseKKTMatrix::Fqq().
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
  /// @brief const version of ImpulseKKTMatrix::Fqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqv() const;

  ///
  /// @brief Jacobian of the state equation of the generalized velocity with 
  /// respect to the impulse forces. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x ImpulseStatus::dimp().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvf();

  ///
  /// @brief const version of ImpulseKKTMatrix::Fvf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvf() const;

  ///
  /// @brief Jacobian of the state equation of the generalized velocity with 
  /// respect to the configuration. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvq();

  ///
  /// @brief const version of ImpulseKKTMatrix::Fvq().
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
  /// @brief const version of ImpulseKKTMatrix::Fvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvv() const;

  ///
  /// @brief Jacobian of the state equation with respect to the impulse forces. 
  /// @return Reference to the block part of the Hessian. 
  /// Size is 2 * Robot::dimv() x ImpulseStatus::dimp().
  ///
  Eigen::Block<Eigen::MatrixXd> Fxf();

  ///
  /// @brief const version of ImpulseKKTMatrix::Fxf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fxf() const;

  ///
  /// @brief Jacobian of the state equation with respect to the state, i.e.,
  /// the configuration and the generalized velocity.
  /// @return Reference to the block part of the Hessian. 
  /// Size is 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fxx();

  ///
  /// @brief const version of ImpulseKKTMatrix::Fxx().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fxx() const;

  ///
  /// @brief Jacobian of the contact position constraint related to impulse 
  /// condition with respect to the configuration. 
  /// ImpulseKKTMatrix::setImpulseStatus() muset be called to set the impulse 
  /// dimension before calling this function.
  /// @return Reference to the block part of the Hessian. 
  /// Size is ImpulseStatus::dimp() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Pq();

  ///
  /// @brief const version of ImpulseKKTMatrix::Pq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Pq() const;

  ///
  /// @brief Jacobian of the contact velocity constraint after impulse with 
  /// respect to the configuration. ImpulseKKTMatrix::setImpulseStatus() must be
  /// called to set the impulse dimension before calling this function.
  /// @return Reference to the block part of the Hessian. 
  /// Size is ImpulseStatus::dimp() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Vq();

  ///
  /// @brief const version of ImpulseKKTMatrix::Vq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Vq() const;

  ///
  /// @brief Jacobian of the contact velocity constraint after impulse with 
  /// respect to the velocity. ImpulseKKTMatrix::setImpulseStatus() must be
  /// called to set the impulse dimension before calling this function.
  /// @return Reference to the block part of the Hessian. 
  /// Size is ImpulseStatus::dimp() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Vv();

  ///
  /// @brief const version of ImpulseKKTMatrix::Vv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Vv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the impulse change in 
  /// velocity and contact forces. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() + ContactStatus::dimf() 
  /// x  Robot::dimv() + ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qdvdvff();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qdvdvff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qdvdvff() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to impulse change in 
  /// velocity. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qdvdv();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qdvdv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qdvdv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration. 
  /// @return Reference to the Hessian. Size is 
  /// ContactStatus::dimf() x ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qff();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qff() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to impulse forces and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is 
  /// ContactStatus::dimf() x ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qfq();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qfq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qfq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to impulse forces and 
  /// velocity. 
  /// @return Reference to the Hessian. Size is 
  /// ContactStatus::dimf() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qfv();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qfv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qfv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and 
  /// impulse forces.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqf();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qqf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqf() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqq();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and 
  /// impulse forces.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvf();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qvf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvf() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity and 
  /// configuration. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to generalized velocity.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxx();

  ///
  /// @brief const version of ImpulseKKTMatrix::Qxx().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxx() const;

  ///
  /// @brief Symmetrize the Hessian for matrix inversion. 
  ///
  void symmetrize();

  ///
  /// @brief Invert the KKT matrix. 
  /// @param[out] KKT_matrix_inverse Inverse of the KKT matrix. Size must 
  /// be ImpulseKKTMatrix::dimKKT() x ImpulseKKTMatrix::dimKKT().
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
  /// @brief Chech the equivalence of two ImpulseKKTMatrix.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const ImpulseKKTMatrix& other) const;

  ///
  /// @brief Chech this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration of the previous time step q_prev.
  Eigen::MatrixXd Fqq_prev;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  DynamicSchurComplement schur_complement_;
  Eigen::MatrixXd FC_, Q_;
  bool has_floating_base_;
  int dimv_, dimx_, dimu_, dim_passive_, dimf_, dimp_, u_begin_, q_begin_, 
      v_begin_, dimKKT_;
  static constexpr int kDimFloatingBase = 6;

};

} // namespace idocp 

#include "idocp/impulse/impulse_kkt_matrix.hxx"

#endif // IDOCP_IMPULSE_KKT_MATRIX_HPP_ 