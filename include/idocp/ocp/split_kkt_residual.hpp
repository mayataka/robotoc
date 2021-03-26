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
  /// @brief Construct a split KKT residual.
  /// @param[in] robot Robot model. 
  ///
  SplitKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~SplitKKTResidual();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitKKTResidual(const SplitKKTResidual&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitKKTResidual& operator=(const SplitKKTResidual&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitKKTResidual(SplitKKTResidual&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
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
  /// @brief Residual in the state equation of q.
  /// @return Reference to the residual in the state equation of q. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fq();

  ///
  /// @brief const version of SplitKKTResidual::Fq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const;

  ///
  /// @brief Residual in the state equation of v.
  /// @return Reference to the residual in the state equation of v. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fv();

  ///
  /// @brief const version of SplitKKTResidual::Fq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const;

  ///
  /// @brief Residual in the state equation.
  /// @return Reference to the residual in the state equation. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fx();

  ///
  /// @brief const version of SplitKKTResidual::Fx().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fx() const;

  ///
  /// @brief Residual in the contact position constraints.
  /// @return Reference to the residual in the contact position constraints. 
  /// Size is ImpulseStatus::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> P();

  ///
  /// @brief const version of SplitKKTResidual::P().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> P() const;

  ///
  /// @brief KKT residual with respect to the control input torques u. 
  /// @return Reference to the residual with respect to u. Size is Robot::dimu().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lu();

  ///
  /// @brief const version of SplitKKTResidual::lu().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lu() const;

  ///
  /// @brief KKT residual with respect to the configuration q. 
  /// @return Reference to the residual with respect to q. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lq();

  ///
  /// @brief const version of SplitKKTResidual::lq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lq() const;

  ///
  /// @brief KKT residual with respect to the velocity v. 
  /// @return Reference to the residual with respect to v. Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lv();

  ///
  /// @brief const version of SplitKKTResidual::lv().
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
  /// @brief Split KKT residual at a time stage. 
  /// @return Reference to the split KKT residual. Size is 
  /// 4 * Robot::dimv() + Robot::dimu().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> splitKKTResidual();

  ///
  /// @brief const version of SplitKKTResidual::splitKKTResidual().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> splitKKTResidual() const;

  ///
  /// @brief KKT residual with respect to the stack of the contact forces f. 
  /// @return Reference to the residual with respect to f. Size is 
  /// SplitKKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief const version of SplitKKTResidual::lf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

  ///
  /// @brief Sets the split KKT residual zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the stack of the contact forces at the 
  /// current contact status.
  /// @return Dimension of the stack of the contact forces.
  ///
  int dimf() const;

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the current 
  /// impulse status.
  /// @return Dimension of the stack of impulse forces.
  ///
  int dimi() const;

  ///
  /// @brief Checks the equivalence of two SplitKKTResidual.
  /// @param[in] other Other object.
  /// @return true if this and other is same. false otherwise.
  ///
  bool isApprox(const SplitKKTResidual& other) const;

  ///
  /// @brief Checks this has at least one NaN.
  /// @return true if this has at least one NaN. false otherwise.
  ///
  bool hasNaN() const;

  /// 
  /// @brief Residual with respect to the acceleration a. Size is Robot::dimv().
  /// 
  Eigen::VectorXd la;

  /// 
  /// @brief Split KKT residual with respect to the virtual control input  
  /// torques on the passive joints.
  /// 
  Vector6d lu_passive;

  /// 
  /// @brief Residual in the part of the state equation.
  /// 
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