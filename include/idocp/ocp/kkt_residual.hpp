#ifndef IDOCP_KKT_RESIDUAL_HPP_
#define IDOCP_KKT_RESIDUAL_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


namespace idocp {

///
/// @class KKTResidual
/// @brief KKT residual split into each time stage. 
///
class KKTResidual {
public:
  ///
  /// @brief Construct a KKT residual.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  KKTResidual(const Robot& robot);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  KKTResidual();

  ///
  /// @brief Destructor. 
  ///
  ~KKTResidual();

  ///
  /// @brief Use default copy constructor. 
  ///
  KKTResidual(const KKTResidual&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  KKTResidual& operator=(const KKTResidual&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  KKTResidual(KKTResidual&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  KKTResidual& operator=(KKTResidual&&) noexcept = default;

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
  /// Size is 2*Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> Fx();

  ///
  /// @brief Residual of the equality constraints.
  /// @return Reference to the residual of the equality constraints. Size is 
  /// KKTResidual::dimc().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> C();

  ///
  /// @brief Residual of the floating base constraint.
  /// @return Reference to the residual of the floating base constraint. Size is 
  /// Robot::dim_passive().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> C_floating_base();

  ///
  /// @brief Residual of the active contact constraints.
  /// @return Reference to the residual of the active contact constraints.  
  /// Size is KKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> C_contacts();

  ///
  /// @brief Residual with respect to acceleration a.
  /// @return Reference to the residual with respect to acceleration a. Size is
  /// Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> la();

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
  /// 2*Robot::dimv().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lx();

  ///
  /// @brief Residual with respect to acceleration and the stack of the 
  /// contact forces, a and f.
  /// @return Reference to the residual with respect to a and f. Size is 
  /// Robot::dimv()+KKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> laf();

  ///
  /// @brief KKT residual.
  /// @return Const reference to the KKT residual. Size is 
  /// KKTResidual::dimKKT().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> KKT_residual() const;

  ///
  /// @brief Residual with respect to q transition.
  /// @return Const reference to the residual with respect to q transition.
  /// Size is Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fq() const;

  ///
  /// @brief Residual with respect to v transition.
  /// @return Const reference to the residual with respect to v transition.
  /// Size is Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fv() const;

  ///
  /// @brief Residual with respect to q and v transition.
  /// @return Const reference to the residual with respect to q and v transition.
  /// Size is 2*Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> Fx() const;

  ///
  /// @brief Residual of the equality constraints.
  /// @return Const reference to the residual of the equality constraints.
  /// Size is KKTResidual::dimc().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> C() const;

  ///
  /// @brief Residual of the floating base constraint.
  /// @return Const reference to the residual of the floating base constraint.  
  /// Size is Robot::dim_passive().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> C_floating_base() const; 

  ///
  /// @brief Residual of the active contact constraints.
  /// @return Const reference to the residual of the active contact constraints.  
  /// Size is KKTResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> C_contacts() const;

  ///
  /// @brief Residual with respect to acceleration a.
  /// @return Const reference to the residual with respect to acceleration a.
  /// Size is Robot::dimv().Size is KKTResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> la() const;

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Const reference to the residual with respect to the stack of the  
  /// contact forces f. Size is KKTResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lf() const;

  ///
  /// @brief Residual with respect to configuration q.
  /// @return Const reference to the residual with respect to configuration q.
  /// Size is Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lq() const;

  ///
  /// @brief Residual with respect to generalized velocity v.
  /// @return Const reference to the residual with respect to generalized 
  /// velocity v. Size is Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lv() const;

  ///
  /// @brief Residual with respect to state x.
  /// @return Const reference to the residual with respect to state x. Size is 
  /// 2*Robot::dimv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> lx() const;

  ///
  /// @brief Residual with respect to acceleration and the stack of the 
  /// contact forces, a and f.
  /// @return Const reference to the residual with respect to a and f. Size is 
  /// Robot::dimv()+KKResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> laf() const;

  ///
  /// @brief Returns the squared norm of the KKT residual.
  /// @param[in] dtau Time step of the horzion.
  ///
  double squaredKKTErrorNorm(const double dtau) const;

  ///
  /// @brief Set the KKT residual zero.
  ///
  void setZeroMinimum();

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

  /// @brief Residual with respect to control input torques u.
  Eigen::VectorXd lu;

  /// @brief Residual of inverse dynamics constraint.
  Eigen::VectorXd u_res;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd kkt_residual_;
  int dimv_, dimx_, dim_passive_, dimf_, dimc_, dimKKT_, max_dimKKT_;

};

} // namespace idocp 

#include "idocp/ocp/kkt_residual.hxx"

#endif // IDOCP_KKT_RESIDUAL_HPP_