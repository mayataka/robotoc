#ifndef IDOCP_IMPULSE_SPLIT_SOLUTION_HPP_ 
#define IDOCP_IMPULSE_SPLIT_SOLUTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

///
/// @class SplitSolution
/// @brief Solution split into each time stage. 
///
class ImpulseSplitSolution {
public:
  ///
  /// @brief Construct a split solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseSplitSolution(const Robot& robot);

  ///
  /// @brief Default constructor. Does not construct any datas. 
  ///
  ImpulseSplitSolution();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitSolution();

  ///
  /// @brief Use default copy constructor. 
  ///
  ImpulseSplitSolution(const ImpulseSplitSolution&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  ImpulseSplitSolution& operator=(const ImpulseSplitSolution&) = default;

  ///
  /// @brief Use default move constructor. 
  ///
  ImpulseSplitSolution(ImpulseSplitSolution&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  ImpulseSplitSolution& operator=(ImpulseSplitSolution&&) noexcept = default;

  ///
  /// @brief Set contact status from robot model, i.e., set dimension of the 
  /// contacts and equality constraints.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Stack of Lagrange multiplier with respect to contact position and 
  /// velocity constraints that is active at the current contact status. Size is 
  /// 2 * ContactStatus::dimf().
  /// @return Reference to the stack of Lagrange multiplier with respect to 
  /// contact position and velocity constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> mu_stack();

  ///
  /// @brief Stack of Lagrange multiplier with respect to contact position 
  /// constraints that is active at the current contact status. Size is 
  /// ContactStatus::dimf().
  /// @return Reference to the stack of Lagrange multiplier with respect to 
  /// contact position constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> mu_stack_contact_position();

  ///
  /// @brief Stack of Lagrange multiplier with respect to contact position and 
  /// velocity constraints that is active at the current contact status. Size is 
  /// 2 * ContactStatus::dimf().
  /// @return Reference to the stack of Lagrange multiplier with respect to 
  /// contact position and velocity constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> mu_stack_contact_velocity();

  ///
  /// @brief Stack of Lagrange multiplier with respect to contact position and 
  /// velocity constraints that is active at the current contact status. Size is 
  /// 2 * ContactStatus::dimf().
  /// @return Const reference to the stack of Lagrange multiplier with respect 
  /// to contact position and velocity constraints.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> mu_stack() const;

  ///
  /// @brief Stack of Lagrange multiplier with respect to contact position 
  /// constraints that is active at the current contact status. Size is 
  /// ContactStatus::dimf().
  /// @return Const reference to the stack of Lagrange multiplier with respect 
  /// to contact position constraints.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> 
  mu_stack_contact_position() const;

  ///
  /// @brief Stack of Lagrange multiplier with respect to contact position and 
  /// velocity constraints that is active at the current contact status. Size is 
  /// 2 * ContactStatus::dimf().
  /// @return Const Reference to the stack of Lagrange multiplier with respect 
  /// to contact position and velocity constraints.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> 
  mu_stack_contact_velocity() const;

  ///
  /// @brief Set the stack of the Lagrange multiplier with respect to active 
  /// equality constraint from mu_floating_base and mu_contacts.
  ///
  void set_mu_stack();

  ///
  /// @brief Set the Lagrange multiplier with respect to active contact 
  /// constraints from mu_stack.
  ///
  void set_mu_contact();

  ///
  /// @brief Stack of active contact forces. Size is Robot::dimf().
  /// @return Reference to the stack of active contact forces.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> f_stack();

  ///
  /// @brief Stack of active contact forces. Size is Robot::dimf().
  /// @return Const reference to the stack of active contact forces.
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> f_stack() const;

  ///
  /// @brief Set the stack of contact forces from each contact forces.
  ///
  void set_f_stack();

  ///
  /// @brief Set the each contact forces from stack of contact forces.
  ///
  void set_f();

  ///
  /// @brief Returns the number of active contacts.
  /// @return Number of active contacts.
  ///
  bool isContactActive(const int contact_index) const;

  ///
  /// @brief Returns the number of active contacts.
  /// @return Number of active contacts.
  ///
  int num_active_contacts() const;

  ///
  /// @brief Returns the dimension of equality constraint at the current 
  /// contact status.
  /// @return Dimension of equality constraint.
  ///
  int dimc() const;

  ///
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of contact forces.
  ///
  int dimf() const;

  ///
  /// @brief Lagrange multiplier with respect to transition of q. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd lmd;

  ///
  /// @brief Lagrange multiplier with respect to transition of v. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd gmm;

  ///
  /// @brief Lagrange multiplier with respect to contact constraint. 
  /// Size is Robot::max_point_contacts().
  ///
  std::vector<Eigen::Vector3d> mu_contact_velocity;

  ///
  /// @brief Lagrange multiplier with respect to contact constraint. 
  /// Size is Robot::max_point_contacts().
  ///
  std::vector<Eigen::Vector3d> mu_contact_position;

  ///
  /// @brief Changes in generalized velocity due to impact. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd dv;

  ///
  /// @brief Contact forces. 
  /// Size is Robot::max_point_contacts().
  ///
  std::vector<Eigen::Vector3d> f;

  ///
  /// @brief Configuration. Size is Robot::dimq().
  ///
  Eigen::VectorXd q;

  ///
  /// @brief Generalized velocity. Size is Robot::dimv().
  ///
  Eigen::VectorXd v;

  ///
  /// @brief Lagrange multiplier with respect to inverse dynamics. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorXd beta;

  ///
  /// @brief Generates split solution filled randomly.
  /// @return Split solution filled randomly.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status.
  ///
  static ImpulseSplitSolution Random(const Robot& robot, 
                                     const ContactStatus& contact_status);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  /// @brief Stack of Lagrange multiplier with respect to equality constraints. 
  Eigen::VectorXd mu_stack_;

  /// @brief Stack of the contact forces. 
  Eigen::VectorXd f_stack_;

  /// @brief Dimension of passive joints. 
  std::vector<bool> is_contact_active_;

  /// @brief Dimension of contact forces at the current contact status. 
  int dimf_;

  /// @brief Dimension of equality constraints at the current contact status. 
  int dimc_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_solution.hxx"

#endif // IDOCP_IMPULSE_SPLIT_SOLUTION_HPP_ 