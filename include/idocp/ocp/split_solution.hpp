#ifndef IDOCP_SPLIT_SOLUTION_HPP_
#define IDOCP_SPLIT_SOLUTION_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

///
/// @class SplitSolution
/// @brief Solution to the optimal control problem split into a time stage. 
///
class SplitSolution {
public:
  ///
  /// @brief Construct a split solution.
  /// @param[in] robot Robot model. 
  ///
  SplitSolution(const Robot& robot);

  ///
  /// @brief Default constructor. 
  ///
  SplitSolution();

  ///
  /// @brief Destructor. 
  ///
  ~SplitSolution();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitSolution(const SplitSolution&) = default;

  ///
  /// @brief Default copy assign operator. 
  ///
  SplitSolution& operator=(const SplitSolution&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  SplitSolution(SplitSolution&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitSolution& operator=(SplitSolution&&) noexcept = default;

  ///
  /// @brief Set contact status, i.e., set the dimension of the contacts.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Set contact status, i.e., set the dimension of the contacts.
  /// @param[in] other Other split solution.
  ///
  void setContactStatus(const SplitSolution& other);

  ///
  /// @brief Set impulse status, i.e., set the dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set impulse status, i.e., set the dimension of the impulse.
  /// @param[in] other Other split solution.
  ///
  void setImpulseStatus(const SplitSolution& other);

  ///
  /// @brief Set impulse status zero.
  ///
  void setImpulseStatus();

  ///
  /// @brief Configuration. Size is Robot::dimq().
  ///
  Eigen::VectorXd q;

  ///
  /// @brief Generalized velocity. Size is Robot::dimv().
  ///
  Eigen::VectorXd v;

  ///
  /// @brief Generalized acceleration. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd a;

  ///
  /// @brief Control input torques. Size is Robot::dimu().
  ///
  Eigen::VectorXd u;

  ///
  /// @brief Contact forces. 
  /// Size is Robot::maxPointContacts().
  ///
  std::vector<Eigen::Vector3d> f;

  ///
  /// @brief Stack of the active contact forces. Size is ContactStatus::dimf().
  /// @return Reference to the stack of the active contact forces.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> f_stack();

  ///
  /// @brief Const version of SplitSolution::f_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> f_stack() const;

  ///
  /// @brief Sets SplitSolution::f_stack() from SplitSolution::f.
  ///
  void set_f_stack();

  ///
  /// @brief Sets SplitSolution::f from SplitSolution::f_stack().
  ///
  void set_f_vector();

  ///
  /// @brief Lagrange multiplier w.r.t. the state equation w.r.t. q.
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd lmd;

  ///
  /// @brief Lagrange multiplier w.r.t. the state equation w.r.t. v.
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd gmm;

  ///
  /// @brief Lagrange multiplier w.r.t. the inverse dynamics constraint. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorXd beta;

  ///
  /// @brief Lagrange multiplier w.r.t. the acceleration-level contact  
  /// constraint. Size is Robot::maxPointContacts().
  ///
  std::vector<Eigen::Vector3d> mu;

  ///
  /// @brief Lagrange multiplier w.r.t. the passive joint constraint. Size is 
  /// Robot::dim_passive().
  ///
  Eigen::VectorXd nu_passive;

  ///
  /// @brief Stack of the Lagrange multipliers w.r.t. the acceleration-level 
  /// contact constraints that is active at the current contact status. Size is 
  /// SplitSolution::dimf().
  /// @return Reference to the stack of the Lagrange multipliers w.r.t.  the 
  /// acceleration-level contact constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> mu_stack();

  ///
  /// @brief Const version of SplitSolution::mu_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> mu_stack() const;

  ///
  /// @brief Set SplitSolution::mu_stack() from SplitSolution::mu. 
  ///
  void set_mu_stack();

  ///
  /// @brief Set SplitSolution::mu from SplitSolution::mu_stack(). 
  ///
  void set_mu_vector();

  ///
  /// @brief Stack of the Lagrange multipliers w.r.t. the switching constraints
  /// that is active at the future impulse status. Size is 
  /// ImpulseSplitSolution::dimf().
  /// @return Reference to the stack of the Lagrange multipliers w.r.t. the
  /// switching constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> xi_stack();

  ///
  /// @brief const version of ImpulseSplitSolution::xi_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> xi_stack() const;

  ///
  /// @brief Returns the dimension of the contact at the current contact status.
  /// @return Dimension of the contact.
  ///
  int dimf() const;

  ///
  /// @brief Returns the dimension of the stack of impulse forces at the 
  /// current impulse status.
  /// @return Dimension of the impulse forces.
  ///
  int dimi() const;

  ///
  /// @brief Return true if there are active contacts and false if not.
  /// @return true if there are active contacts and false if not. 
  ///
  bool hasActiveContacts() const;

  ///
  /// @brief Returns true if there are active impulse constraints and false if 
  /// not.
  /// @return true if there are active impulse constraints and false if not. 
  ///
  bool hasActiveImpulse() const;

  ///
  /// @brief Return true if contact is active and false if not.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return true if a contact is active and false if not. 
  ///
  bool isContactActive(const int contact_index) const;

  ///
  /// @brief Return activities of contacts.
  /// @return Activities of contacts. 
  ///
  std::vector<bool> isContactActive() const;

  ///
  /// @brief Integrates the solution based on step size and direction. 
  /// @param[in] robot Robot model.
  /// @param[in] step_size Step size.
  /// @param[in] d Split direction.
  ///
  void integrate(const Robot& robot, const double step_size, 
                 const SplitDirection& d);

  ///
  /// @brief Copies the primal solution from another solution. 
  /// @param[in] other Another split solution.
  ///
  void copyPrimal(const SplitSolution& other);

  ///
  /// @brief Copies the dual solution from another solution. 
  /// @param[in] other Another split solution.
  ///
  void copyDual(const SplitSolution& other);

  ///
  /// @brief Return true if two SplitSolution have the same value and false if 
  /// not. 
  /// @param[in] other Split solution that is compared with this object.
  ///
  bool isApprox(const SplitSolution& other) const;

  ///
  /// @brief Set each component vector by random value based on the current 
  /// contact status. 
  /// @param[in] robot Robot model.
  ///
  void setRandom(const Robot& robot);

  ///
  /// @brief Set each component vector by random value. Contact status is reset.
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  ///
  void setRandom(const Robot& robot, const ContactStatus& contact_status);

  ///
  /// @brief Set each component vector by random value. Impulse status is reset.
  /// @param[in] robot Robot model.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const Robot& robot, const ImpulseStatus& impulse_status);

  ///
  /// @brief Set each component vector by random value. Contact status and 
  /// impulse status are reset.
  /// @param[in] robot Robot model.
  /// @param[in] contact_status Contact status.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const Robot& robot, const ContactStatus& contact_status,
                 const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates split solution filled randomly.
  /// @return Split solution filled randomly.
  /// @param[in] robot Robot model. 
  ///
  static SplitSolution Random(const Robot& robot);

  ///
  /// @brief Generates split solution filled randomly.
  /// @return Split solution filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status.
  ///
  static SplitSolution Random(const Robot& robot, 
                              const ContactStatus& contact_status);

  ///
  /// @brief Generates split solution filled randomly.
  /// @return Split solution filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status.
  ///
  static SplitSolution Random(const Robot& robot, 
                              const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates split solution filled randomly.
  /// @return Split solution filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] contact_status Contact status.
  /// @param[in] impulse_status Impulse status.
  ///
  static SplitSolution Random(const Robot& robot, 
                              const ContactStatus& contact_status,
                              const ImpulseStatus& impulse_status);

  ///
  /// @brief Displays the split solution onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const SplitSolution& s);

private:
  Eigen::VectorXd mu_stack_, f_stack_, xi_stack_;
  bool has_floating_base_, has_active_contacts_, has_active_impulse_;
  std::vector<bool> is_contact_active_;
  int dimf_, dimi_;

};

} // namespace idocp 

#include "idocp/ocp/split_solution.hxx"

#endif // IDOCP_SPLIT_SOLUTION_HPP_