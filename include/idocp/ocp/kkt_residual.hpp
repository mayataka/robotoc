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
  /// @brief Residual with respect to acceleration and the stack of the 
  /// contact forces, a and f.
  /// @return Reference to the residual with respect to a and f. Size is 
  /// Robot::dimv() + KKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lu();

  ///
  /// @brief Residual with respect to acceleration and the stack of the 
  /// contact forces, a and f.
  /// @return Reference to the residual with respect to a and f. Size is 
  /// Robot::dimv() + KKTResidual::dimf().
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

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is KKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> C();

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is KKTResidual::dimf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> C() const;

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is KKTResidual::dimf().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> lf();

  ///
  /// @brief Residual with respect to the stack of the contact forces f.
  /// @return Reference to the residual with respect to the stack of the  
  /// contact forces f. Size is KKTResidual::dimf().
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

  /// @brief KKT residual.
  Eigen::VectorXd KKT_residual;

  /// @brief Residual with respect to control input torques u.
  Eigen::VectorXd la;

  /// @brief Residual of inverse dynamics constraint.
  Eigen::VectorXd ID;

  /// @brief Residual with respect to control input torques u.
  Eigen::VectorXd lu_passive;

  /// @brief Residual of floating base constraint.
  Eigen::VectorXd C_passive;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd C_full_, lf_full_;
  int dimv_, dimx_, dimu_, dim_passive_, dimf_, dimKKT_;

};

} // namespace idocp 

#include "idocp/ocp/kkt_residual.hxx"

#endif // IDOCP_KKT_RESIDUAL_HPP_