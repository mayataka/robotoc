#ifndef IDOCP_KKT_MATRIX_HPP_
#define IDOCP_KKT_MATRIX_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
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
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqu();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqq();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqv();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvu();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvq();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvv();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fxu();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fxu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fxx();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fxx() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_full();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_passive_topLeft();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_passive_topLeft() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_passive_topRight();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_passive_topRight() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu_passive_bottomLeft();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu_passive_bottomLeft() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quu();

  ///
  /// @brief Hessian of the Lagrangian with respect to the control input torques 
  /// u. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq_full();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq_passive();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quq();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv_full();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv_passive();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// configuration, a and q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Quv();

  ///
  /// @brief Hessian of the Lagrangian with respect to acceleration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Quv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu_full();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu_passive();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqq();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief Hessian of the Lagrangian with respect to configuration and 
  /// velocity, a and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu_full();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu_passive();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvq() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief Hessian of the Lagrangian with respect to velocity and v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvv() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qux_full();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qux_full() const;
  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qux_passive();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qux_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qux();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qux() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxu_full();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxu_full() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxu_passive();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxu_passive() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxu();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qxx();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qxx() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaaff();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaaff() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qaa();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qaa() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qff();

  ///
  /// @brief Hessian of the Lagrangian with respect to state, q and v. 
  /// @return Reference to the Hessian. Size is 
  /// 2 * Robot::dimv() x 2 * Robot::dimv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qff() const;

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
  void invert(const Eigen::MatrixBase<MatrixType>& KKT_matrix_inverse);

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
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of the stack of contact forces.
  ///
  int dimf() const;

  bool isApprox(const KKTMatrix& other) const;

  /// @brief Derivative of the state equation with respect to the 
  /// configuration of the previous time step q_prev.
  Eigen::MatrixXd Fqq_prev;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  SchurComplement schur_complement_;
  Eigen::MatrixXd F_, Q_, Qaaff_full_;
  bool has_floating_base_;
  int dimv_, dimx_, dimu_, dim_passive_, dimf_, u_begin_, q_begin_, v_begin_, 
      dimQ_, dimKKT_;
  static constexpr int kDimFloatingBase = 6;

};

} // namespace idocp 

#include "idocp/ocp/kkt_matrix.hxx"

#endif // IDOCP_KKT_MATRIX_HPP_