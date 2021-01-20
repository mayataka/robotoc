#ifndef IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_ 
#define IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class ImpulseSplitKKTResidual
/// @brief KKT residual split at the impulse stage. 
///
class ImpulseSplitKKTResidual {
public:
  ///
  /// @brief Construct a KKT residual.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseSplitKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  ImpulseSplitKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitKKTResidual();

  ///
  /// @brief Use default copy constructor. 
  ///
  ImpulseSplitKKTResidual(const ImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  ImpulseSplitKKTResidual& operator=(const ImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  ImpulseSplitKKTResidual(ImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  ImpulseSplitKKTResidual& operator=(
      ImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Set impulse status from robot model, i.e., set dimension of the 
  /// impulse and equality constraints.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

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
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> V();

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> V() const;

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is ImpulseSplitKKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is ImpulseSplitKKTResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

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

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> P();

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> P() const;

  Eigen::VectorBlock<Eigen::VectorXd> splitKKTResidual();

  const Eigen::VectorBlock<const Eigen::VectorXd> splitKKTResidual() const;

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
  /// @brief Chech the equivalence of two ImpulseSplitKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const ImpulseSplitKKTResidual& other) const;

  ///
  /// @brief Chech this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  /// @brief Residual with respect to control input torques u.
  Eigen::VectorXd ldv;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd KKT_residual_full_, P_full_;
  int dimv_, dimx_, dimf_, dimKKT_, lf_begin_, lq_begin_, lv_begin_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_kkt_residual.hxx"

#endif // IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_ 