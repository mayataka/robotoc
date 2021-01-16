#ifndef IDOCP_PRE_IMPULSE_SPLIT_KKT_MATRIX_HPP_
#define IDOCP_PRE_IMPULSE_SPLIT_KKT_MATRIX_HPP_

#include "Eigen/Core"
#include "Eigen/LU"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class PreImpulseSplitKKTMatrix
/// @brief The KKT matrix split into a time stage.
///
class PreImpulseSplitKKTMatrix {
public:
  ///
  /// @brief Construct a KKT matrix.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  PreImpulseSplitKKTMatrix(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  PreImpulseSplitKKTMatrix();

  ///
  /// @brief Destructor. 
  ///
  ~PreImpulseSplitKKTMatrix();

  ///
  /// @brief Default copy constructor. 
  ///
  PreImpulseSplitKKTMatrix(const PreImpulseSplitKKTMatrix&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  PreImpulseSplitKKTMatrix& operator=(const PreImpulseSplitKKTMatrix&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  PreImpulseSplitKKTMatrix(PreImpulseSplitKKTMatrix&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  PreImpulseSplitKKTMatrix& operator=(PreImpulseSplitKKTMatrix&&) noexcept = default;

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] impulse_status Contact status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

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
  /// @brief const version of ImpulseSplitKKTMatrix::Pq().
  ///
  const Eigen::Block<const Eigen::MatrixXd> Pq() const;

  ///
  /// @brief Set the all components zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of the stack of impulse forces.
  ///
  int dimf() const;

  ///
  /// @brief Chech the equivalence of two PreImpulseSplitKKTMatrix.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const PreImpulseSplitKKTMatrix& other) const;

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
  Eigen::MatrixXd Pq_full_;
  int dimv_, dimf_;

};

} // namespace idocp 

#include "idocp/ocp/pre_impulse_split_kkt_matrix.hxx"

#endif // IDOCP_PRE_IMPULSE_SPLIT_KKT_MATRIX_HPP_ 