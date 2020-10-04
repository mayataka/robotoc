#ifndef IDOCP_SPLIT_DIRECTION_HPP_
#define IDOCP_SPLIT_DIRECTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitDirection
/// @brief Newton direction split into each time stage. 
///
class SplitDirection {
public:
  ///
  /// @brief Construct a split solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  SplitDirection(const Robot& robot);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  SplitDirection();

  ///
  /// @brief Destructor. 
  ///
  ~SplitDirection();

  ///
  /// @brief Use default copy constructor. 
  ///
  SplitDirection(const SplitDirection&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  SplitDirection& operator=(const SplitDirection&) = default;
 
  ///
  /// @brief Use default move constructor. 
  ///
  SplitDirection(SplitDirection&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  SplitDirection& operator=(SplitDirection&&) noexcept = default;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Returns the stack of Newton directions. Size is 
  /// SplitDirection::dimKKT().
  /// @return Reference to the stack of Newton directions.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> split_direction();

  ///
  /// @brief Returns the Newton direction of lmd. Size is Robot::dimv().
  /// @return Reference to the Newton direction of lmd.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dlmd();

  ///
  /// @brief Returns the Newton direction of gmm. Size is Robot::dimv().
  /// @return Reference to the Newton direction of gmm.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dgmm();

  ///
  /// @brief Returns the Newton direction of mu_stack. Size is 
  /// Robot::dim_passive() + Robot::dimf().
  /// @return Reference to the Newton direction of mu_stack.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dmu();

  ///
  /// @brief Returns the Newton direction of a. Size is Robot::dimv().
  /// @return Reference to the Newton direction of a.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> da();

  ///
  /// @brief Returns the Newton direction of f_stack. Size is Robot::dimf().
  /// @return Reference to the Newton direction of f_stack.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> df();

  ///
  /// @brief Returns the Newton direction of q. Size is Robot::dimv().
  /// @return Reference to the Newton direction of q.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dq();

  ///
  /// @brief Returns the Newton direction of v. Size is Robot::dimv().
  /// @return Reference to the Newton direction of v.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dv();

  ///
  /// @brief Returns the Newton direction of q and v. Size is 2 * Robot::dimv().
  /// @return Reference to the Newton direction of q and v.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dx();

  ///
  /// @brief Returns the stack of Newton directions. Size is 
  /// SplitDirection::dimKKT().
  /// @return Const reference to the stack of Newton directions.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> split_direction() const;

  ///
  /// @brief Returns the Newton direction of lmd. Size is Robot::dimv().
  /// @return Const reference to the Newton direction of lmd.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dlmd() const;

  ///
  /// @brief Returns the Newton direction of gmm. Size is Robot::dimv().
  /// @return Const reference to the Newton direction of gmm.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dgmm() const;

  ///
  /// @brief Returns the Newton direction of mu_stack. Size is 
  /// Robot::dim_passive() + Robot::dimf().
  /// @return Const reference to the Newton direction of mu_stack.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dmu() const;

  ///
  /// @brief Returns the Newton direction of a. Size is Robot::dimv().
  /// @return Const reference to the Newton direction of a.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> da() const;

  ///
  /// @brief Returns the Newton direction of f_stack. Size is Robot::dimf().
  /// @return Const reference to the Newton direction of f_stack.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> df() const;

  ///
  /// @brief Returns the Newton direction of q. Size is Robot::dimv().
  /// @return Const reference to the Newton direction of q.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dq() const;

  ///
  /// @brief Returns the Newton direction of v. Size is Robot::dimv().
  /// @return Const reference to the Newton direction of v.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dv() const;

  ///
  /// @brief Returns the Newton direction of q and v. Size is 2 * Robot::dimv().
  /// @return Const reference to the Newton direction of q and v.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dx() const;

  ///
  /// @brief Set the all alements of the direction to zero.
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

  /// @brief Newton direction of u.
  Eigen::VectorXd du;

  /// @brief Newton direction of beta.
  Eigen::VectorXd dbeta;

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  ///
  static SplitDirection Random(const Robot& robot);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  ///
  static SplitDirection Random(const Robot& robot, 
                               const ContactStatus& contact_status);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  /// @brief Stack of the Newton directions.
  Eigen::VectorXd split_direction_;

  /// @brief Dimension of velocity v.
  int dimv_;

  /// @brief Dimension of state x.
  int dimx_;

  /// @brief Dimension of passive joints. 
  int dim_passive_;

  /// @brief Dimension of contact forces at the current contact status. 
  int dimf_;

  /// @brief Dimension of equality constraints at the current contact status. 
  int dimc_;

  /// @brief Dimension of the split KKT system. 
  int dimKKT_;

  /// @brief Maximum dimension of the split KKT system. 
  int max_dimKKT_;

};

} // namespace idocp 

#include "idocp/ocp/split_direction.hxx"

#endif // IDOCP_SPLIT_OCP_DIRECTION_HPP_