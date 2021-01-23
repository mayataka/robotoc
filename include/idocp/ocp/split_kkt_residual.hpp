#ifndef IDOCP_SPLIT_KKT_RESIDUAL_HPP_ 
#define IDOCP_SPLIT_KKT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class SplitKKTResidual
/// @brief KKT residual split into each time stage. 
///
class SplitKKTResidual {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Construct a KKT residual.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  SplitKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  SplitKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~SplitKKTResidual();

  ///
  /// @brief Use default copy constructor. 
  ///
  SplitKKTResidual(const SplitKKTResidual&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  SplitKKTResidual& operator=(const SplitKKTResidual&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  SplitKKTResidual(SplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  SplitKKTResidual& operator=(SplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Set contact status, i.e., set dimension of the contact.
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
  /// @brief Residual with respect to q transition.
  /// @return Reference to the residual with respect to q transition. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fq();

  ///
  /// @brief Residual with respect to q transition.
  /// @return Reference to the residual with respect to q transition. Size is 
  /// Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const;

  ///
  /// @brief Residual with respect to v transition.
  /// @return Reference to the residual with respect to v transition. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fv();

  ///
  /// @brief Residual with respect to v transition.
  /// @return Reference to the residual with respect to v transition. Size is 
  /// Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const;

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fx();

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fx() const;

  ///
  /// @brief Residual with respect to impulse condition constraint.
  /// @return Reference to the impulse condition constraint.
  /// Size is ImpulseStatus::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> P();

  ///
  /// @brief const version of SplitKKTResidual::P().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> P() const;

  ///
  /// @brief Residual with respect to acceleration and the stack of the 
  /// contact forces, a and f.
  /// @return Reference to the residual with respect to a and f. Size is 
  /// Robot::dimv() + SplitKKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lu();

  ///
  /// @brief Residual with respect to acceleration and the stack of the 
  /// contact forces, a and f.
  /// @return Reference to the residual with respect to a and f. Size is 
  /// Robot::dimv() + SplitKKTResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lu() const;

  ///
  /// @brief Residual with respect to configuration q.
  /// @return Reference to the residual with respect to configuration q. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lq();

  ///
  /// @brief Residual with respect to configuration q.
  /// @return Reference to the residual with respect to configuration q. Size is 
  /// Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lq() const;

  ///
  /// @brief Residual with respect to generalized velocity v.
  /// @return Reference to the residual with respect to generalized velocity v.
  /// Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lv();

  ///
  /// @brief Residual with respect to generalized velocity v.
  /// @return Reference to the residual with respect to generalized velocity v.
  /// Size is Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lv() const;

  ///
  /// @brief Residual with respect to state x.
  /// @return Reference to the residual with respect to state x. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lx();

  ///
  /// @brief Residual with respect to state x.
  /// @return Reference to the residual with respect to state x. Size is 
  /// 2 * Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lx() const;

  Eigen::VectorBlock<Eigen::VectorXd> splitKKTResidual();

  const Eigen::VectorBlock<const Eigen::VectorXd> splitKKTResidual() const;

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is SplitKKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is SplitKKTResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

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

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of the stack of impulse forces.
  ///
  int dimi() const;

  ///
  /// @brief Chech the equivalence of two SplitKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitKKTResidual& other) const;

  ///
  /// @brief Chech this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  /// @brief Residual with respect to control input torques u.
  Eigen::VectorXd la;

  /// @brief Residual with respect to control input torques u.
  Vector6d lu_passive;

  Vector6d Fq_prev;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd kkt_residual_full_, lf_full_;
  int dimv_, dimx_, dimu_, dim_passive_, dimf_, dimi_, dimKKT_,
      lu_begin_, lq_begin_, lv_begin_;
  bool has_floating_base_;

};

} // namespace idocp 

#include "idocp/ocp/split_kkt_residual.hxx"

#endif // IDOCP_SPLIT_KKT_RESIDUAL_HPP_