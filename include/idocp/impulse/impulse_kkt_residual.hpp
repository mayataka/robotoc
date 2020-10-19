#ifndef IDOCP_IMPULSE_KKT_RESIDUAL_HPP_
#define IDOCP_IMPULSE_KKT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


namespace idocp {

///
/// @class ImpulseKKTResidual
/// @brief KKT residual split into impulse time stage. 
///
class ImpulseKKTResidual {
public:
  ///
  /// @brief Construct a impulse KKT residual.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] use_contact_position_constraint true if you treat the contact 
  /// position constraint in this impulse stage. false if not.
  ///
  ImpulseKKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  ImpulseKKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseKKTResidual();

  ///
  /// @brief Use default copy constructor. 
  ///
  ImpulseKKTResidual(const ImpulseKKTResidual&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  ImpulseKKTResidual& operator=(const ImpulseKKTResidual&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  ImpulseKKTResidual(ImpulseKKTResidual&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  ImpulseKKTResidual& operator=(ImpulseKKTResidual&&) noexcept = default;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief KKT residual. 
  /// @return Reference to the KKT residual. Size is KKTResidual::dimKKT().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> KKT_residual();

  ///
  /// @brief Residual with respect to q transition.
  /// @return Reference to the residual with respect to q transition. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fq();

  ///
  /// @brief Residual with respect to v transition.
  /// @return Reference to the residual with respect to v transition. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fv();

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fx();

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Const reference to the residual with respect to q and v transition.
  /// Size is 2 * Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fx() const;

  ///
  /// @brief Residual of the contact position and velocity constraint.
  /// @return Reference to the residual of the contact position constraint.  
  /// Size is 2 * KKTResidual::dimf(). 
  ///
  Eigen::VectorBlock<Eigen::VectorXd> C();

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is KKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief Residual with respect to configuration q.
  /// @return Reference to the residual with respect to configuration q. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lq();

  ///
  /// @brief Residual with respect to generalized velocity v.
  /// @return Reference to the residual with respect to generalized velocity v.
  /// Size is Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lv();

  ///
  /// @brief Residual with respect to state x.
  /// @return Reference to the residual with respect to state x. Size is 
  /// 2 * Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lx();

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
  /// @brief Returns the maximum dimension of the KKT.
  /// @return Maximum dimension of the KKT at the current contact status.
  ///
  int max_dimKKT() const;

  ///
  /// @brief Returns the dimension of equality constraint at the current 
  /// contact status.
  /// @return Dimension of equality constraint.
  ///
  int dimc() const;

  ///
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of the stack of contact forces.
  ///
  int dimf() const;

  /// @brief Residual with respect to the impulse velocity.
  Eigen::VectorXd ldv;

  /// @brief Residual of inverse impulse dynamics constraint.
  Eigen::VectorXd dv_res;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd kkt_residual_;
  int dimv_, dimx_, dimf_, dimc_, dimKKT_, max_dimKKT_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_kkt_residual.hxx"

#endif // IDOCP_IMPULSE_KKT_RESIDUAL_HPP_ 