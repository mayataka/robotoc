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
  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Construct a KKT residual.
  /// @param[in] robot Robot model. 
  ///
  ImpulseSplitKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseSplitKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitKKTResidual();

  ///
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitKKTResidual(const ImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitKKTResidual& operator=(const ImpulseSplitKKTResidual&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitKKTResidual(ImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitKKTResidual& operator=(
      ImpulseSplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Residual in the state equation of q.
  /// @return Reference to the residual in the state equation of q. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fq();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::Fq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const;

  ///
  /// @brief Residual in the state equation of v.
  /// @return Reference to the residual in the state equation of v. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fv();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::Fv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const;

  ///
  /// @brief Residual in the state equation.
  /// @return Reference to the residual in the state equation. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fx();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::Fx().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fx() const;

  ///
  /// @brief Residual in the impulse velocity constraints.
  /// @return Reference to the residual in the impulse velocity constraints. 
  /// Size is ImpulseStatus::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> V();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::V().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> V() const;

  ///
  /// @brief KKT residual with respect to the stack of the impulse forces f. 
  /// @return Reference to the residual with respect to f. Size is 
  /// ImpulseSplitKKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::lf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

  ///
  /// @brief KKT residual with respect to the configuration q. 
  /// @return Reference to the residual with respect to q. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lq();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::lq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lq() const;

  ///
  /// @brief KKT residual with respect to the velocity v. 
  /// @return Reference to the residual with respect to v. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lv();

  ///
  /// @brief const version of ImpulseSplitKKTResidual::lv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lv() const;

  ///
  /// @brief KKT residual with respect to the state x. 
  /// @return Reference to the residual with respect to x. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lx();

  ///
  /// @brief const version of SplitKKTResidual::lx().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lx() const;

  ///
  /// @brief Split KKT residual at a impulse time stage. 
  /// @return Reference to the split KKT residual. Size is 
  /// 4 * Robot::dimv() + ImpulseStatus::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> splitKKTResidual();

  ///
  /// @brief const version of SplitKKTResidual::splitKKTResidual().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> splitKKTResidual() const;

  ///
  /// @brief Sets the split KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of the impluse forces at the 
  /// current impulse status.
  /// @return Dimension of the stack of the impulse forces.
  ///
  int dimf() const;

  ///
  /// @brief Checks the equivalence of two ImpulseSplitKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const ImpulseSplitKKTResidual& other) const;

  ///
  /// @brief Checks his has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  /// 
  /// @brief Residual with respect to the impulse change in the velocity dv. 
  /// Size is Robot::dimv().
  /// 
  Eigen::VectorXd ldv;

  /// 
  /// @brief Residual in the part of the state equation.
  /// 
  Vector6d Fq_prev;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd KKT_residual_full_;
  int dimv_, dimx_, dimf_, dimKKT_, lf_begin_, lq_begin_, lv_begin_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_kkt_residual.hxx"

#endif // IDOCP_IMPULSE_SPLIT_KKT_RESIDUAL_HPP_ 