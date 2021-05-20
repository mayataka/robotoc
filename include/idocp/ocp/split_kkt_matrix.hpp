#ifndef IDOCP_SPLIT_KKT_MATRIX_HPP_
#define IDOCP_SPLIT_KKT_MATRIX_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


namespace idocp {

///
/// @class SplitKKTMatrix
/// @brief The KKT matrix split into a time stage.
///
class SplitKKTMatrix {
public:
  ///
  /// @brief Construct a split KKT matrix.
  /// @param[in] robot Robot model. 
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
  /// @brief Jacobian of the state equation w.r.t. the state x.
  ///
  Eigen::MatrixXd Fxx;

  ///
  /// @brief Jacobian of the state equation (w.r.t. q) w.r.t. q.
  /// @return Reference to the block of the Jacobian of the constraints. Size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqq();

  ///
  /// @brief const version of SplitKKTMatrix::Fqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqq() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. q) w.r.t. v.
  /// @return Reference to the block of the Jacobian of the constraints. Size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqv();

  ///
  /// @brief const version of SplitKKTMatrix::Fqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqv() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. v) w.r.t. q.
  /// @return Reference to the block of the Jacobian of the constraints. Size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvq();

  ///
  /// @brief const version of SplitKKTMatrix::Fvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvq() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. v) to v.
  /// @return Reference to the block of the Jacobian of the constraints. Size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvv();

  ///
  /// @brief const version of SplitKKTMatrix::Fvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvv() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. v) w.r.t. u. 
  ///
  Eigen::MatrixXd Fvu;

  ///
  /// @brief Hessian w.r.t. to the state x and state x.
  ///
  Eigen::MatrixXd Qxx;

  ///
  /// @brief Hessian w.r.t. the configuration q and configuration q.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqq();

  ///
  /// @brief const version of SplitKKTMatrix::Qqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqq() const;

  ///
  /// @brief Hessian w.r.t. the configuration q and joint velocity v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief const version of SplitKKTMatrix::Qqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqv() const;

  ///
  /// @brief Hessian w.r.t. the joint velocity v and configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief const version of SplitKKTMatrix::Qvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvq() const;

  ///
  /// @brief Hessian w.r.t. the joint velocity v and joint velocity v.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief const version of SplitKKTMatrix::Qvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvv() const;

  ///
  /// @brief Hessian w.r.t. the acceleration a and acceleration a. 
  ///
  Eigen::MatrixXd Qaa;

  ///
  /// @brief Hessian w.r.t. the state x and the control input torques u.
  ///
  Eigen::MatrixXd Qxu;

  Eigen::MatrixXd Qxu_passive;

  ///
  /// @brief Hessian of the Lagrangian with respect to the configuration q and
  /// control input torques u. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqu();

  ///
  /// @brief const version of SplitKKTMatrix::Qqu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqu() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the velocity v and
  /// control input torques u. 
  /// @return Reference to the Hessian. Size is Robot::dimu() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvu();

  ///
  /// @brief const version of SplitKKTMatrix::Qvu().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvu() const;

  ///
  /// @brief Hessian w.r.t. the control input torques u and the control input 
  /// torques u.
  ///
  Eigen::MatrixXd Quu;

  Eigen::MatrixXd Quu_passive_topRight;

  ///
  /// @brief Hessian of the Lagrangian with respect to the contact forces f. 
  /// @return Reference to the Hessian. Size is 
  /// ContactStatus::dimf() x ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qff();

  ///
  /// @brief const version of SplitKKTMatrix::Qff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qff() const;

  ///
  /// @brief Hessian of the Lagrangian with respect to the configuration and 
  /// contact forces. 
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x ContactStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqf();

  ///
  /// @brief const version of SplitKKTMatrix::Qqf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqf() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. q) w.r.t. q_prev.
  /// If Robot::hasFloatingBase() is true, size is Robot::dimv() x Robot::dimv().
  /// Otherwise, 0 x 0.
  ///
  Eigen::MatrixXd Fqq_prev;

  ///
  /// @brief Inverse of the Jacobian of the state equation (w.r.t. q) w.r.t. q.
  /// If Robot::hasFloatingBase() is true, size is 6 x 6. Otherwise, 0 x 0.
  ///
  Eigen::MatrixXd Fqq_inv;

  ///
  /// @brief Inverse of the Jacobian of the state equation (w.r.t. q) w.r.t. 
  /// q_prev. If Robot::hasFloatingBase() is true, size is 6 x 6. 
  /// Otherwise, 0 x 0.
  ///
  Eigen::MatrixXd Fqq_prev_inv;

  ///
  /// @brief Set the all components zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of the stack of contact forces.
  ///
  int dimf() const;

  ///
  /// @brief Checks dimensional consistency of each component. 
  /// @return true if the dimension is consistent. false if not.
  ///
  bool isDimensionConsistent() const;

  ///
  /// @brief Checks the equivalence of two SplitKKTMatrix.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitKKTMatrix& other) const;

  ///
  /// @brief Checks this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

private:
  Eigen::MatrixXd Qff_full_, Qqf_full_;
  bool has_floating_base_;
  int dimv_, dimx_, dimu_, dim_passive_, dimf_;

};

} // namespace idocp 

#include "idocp/ocp/split_kkt_matrix.hxx"

#endif // IDOCP_SPLIT_KKT_MATRIX_HPP_ 