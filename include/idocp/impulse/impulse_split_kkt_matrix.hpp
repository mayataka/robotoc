#ifndef IDOCP_IMPULSE_SPLIT_KKT_MATRIX_HPP_
#define IDOCP_IMPULSE_SPLIT_KKT_MATRIX_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class ImpulseSplitKKTMatrix
/// @brief KKT matrix split into an impulse time stage.
///
class ImpulseSplitKKTMatrix {
public:
  ///
  /// @brief Construct a split KKT matrix.
  /// @param[in] robot Robot model. 
  ///
  ImpulseSplitKKTMatrix(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitKKTMatrix();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitKKTMatrix();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitKKTMatrix(const ImpulseSplitKKTMatrix&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitKKTMatrix& operator=(const ImpulseSplitKKTMatrix&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitKKTMatrix(ImpulseSplitKKTMatrix&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitKKTMatrix& operator=(ImpulseSplitKKTMatrix&&) noexcept = default;

  ///
  /// @brief Set the impulse status, i.e., set dimension of the impulses.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

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
  /// @brief const version of ImpulseSplitKKTMatrix::Fqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqq() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. q) w.r.t. v.
  /// @return Reference to the block of the Jacobian of the constraints. Size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fqv();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Fqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fqv() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. v) w.r.t. q.
  /// @return Reference to the block of the Jacobian of the constraints. Size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvq();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Fvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvq() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. v) to v.
  /// @return Reference to the block of the Jacobian of the constraints. Size 
  /// is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Fvv();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Fvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Fvv() const;

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
  /// @brief const version of ImpulseSplitKKTMatrix::Qqq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqq() const;

  ///
  /// @brief Hessian w.r.t. the configuration q and joint velocity v. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqv();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Qqv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqv() const;

  ///
  /// @brief Hessian w.r.t. the joint velocity v and configuration q. 
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvq();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Qvq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvq() const;

  ///
  /// @brief Hessian w.r.t. the joint velocity v and joint velocity v.
  /// @return Reference to the Hessian. Size is Robot::dimv() x Robot::dimv().
  ///
  Eigen::Block<Eigen::MatrixXd> Qvv();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Qvv().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qvv() const;

  ///
  /// @brief Hessian w.r.t. the impulse change in the velocity ddv. 
  ///
  Eigen::MatrixXd Qdvdv;

  ///
  /// @brief Hessian w.r.t. the impulse forces f and impulse forces f.
  /// @return Reference to the Hessian. Size is 
  /// ImpulseStatus::dimf() x ImpulseStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qff();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Qff().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qff() const;

  ///
  /// @brief Hessian w.r.t. the configuration q and impulse forces f.
  /// @return Reference to the Hessian. Size is 
  /// Robot::dimv() x ImpulseStatus::dimf().
  ///
  Eigen::Block<Eigen::MatrixXd> Qqf();

  ///
  /// @brief const version of ImpulseSplitKKTMatrix::Qqf().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Qqf() const;

  ///
  /// @brief Jacobian of the state equation (w.r.t. q) w.r.t. q_prev.
  /// If Robot::hasFloatingBase() is true, size is Robot::dimv() x Robot::dimv().
  /// Otherwise, 0 x 0.
  ///
  Eigen::MatrixXd Fqq_prev;

  ///
  /// @brief Set the all components zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of the stack of impulse forces.
  ///
  int dimi() const;

  ///
  /// @brief Checks dimensional consistency of each component. 
  /// @return true if the dimension is consistent. false if not.
  ///
  bool isDimensionConsistent() const;

  ///
  /// @brief Checks the equivalence of two ImpulseSplitKKTMatrix.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const ImpulseSplitKKTMatrix& other) const;

  ///
  /// @brief Checks if this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  ///
  /// @brief Set by random value based on the current impulse status. 
  ///
  void setRandom();

  ///
  /// @brief Set by random value. Impulse status is reset.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates impulse split KKT matrix filled randomly.
  /// @return Impulse split KKT matrix filled randomly.
  /// @param[in] robot Robot model. 
  ///
  static ImpulseSplitKKTMatrix Random(const Robot& robot);

  ///
  /// @brief Generates impulse split KKT matrix filled randomly.
  /// @return Impulse split KKT matrix filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status.
  ///
  static ImpulseSplitKKTMatrix Random(const Robot& robot, 
                                      const ImpulseStatus& impulse_status);

private:
  Eigen::MatrixXd Qff_full_, Qqf_full_;
  int dimv_, dimi_;
  bool has_floating_base_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_kkt_matrix.hxx"

#endif // IDOCP_IMPULSE_SPLIT_KKT_MATRIX_HPP_ 