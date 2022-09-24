#ifndef ROBOTOC_SPLIT_SOLUTION_HPP_
#define ROBOTOC_SPLIT_SOLUTION_HPP_

#include <vector>
#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/core/split_direction.hpp"


namespace robotoc {

///
/// @class SplitSolution
/// @brief Solution to the optimal control problem split into a time stage. 
///
class SplitSolution {
public:
  using Vector6d = Eigen::Matrix<double, 6, 1>;

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
  /// @brief Default destructor. 
  ///
  ~SplitSolution() = default;

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
  /// @brief Set contact status, i.e., set the dimension of the contacts.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ImpulseStatus& contact_status);

  ///
  /// @brief Sets the dimension of the switching constraint.
  /// @param[in] dims The dimension of the switching constraint. Must be non-negative.
  ///
  void setSwitchingConstraintDimension(const int dims);

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
  /// @brief Generalized acceleration. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd dv;

  ///
  /// @brief Control input torques. Size is Robot::dimu().
  ///
  Eigen::VectorXd u;

  ///
  /// @brief Contact wrenches. Upper 3 elements are linear contact force
  /// and the lower 3 elements are the angular momentum.
  /// Size is Robot::maxNumContacts().
  ///
  std::vector<Vector6d> f;

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
  /// constraint. Upper 3 elements are w.r.t. the linear contact acceleration
  /// and the lower 3 elements are w.r.t. the angular contact acceleration.
  /// Size is Robot::maxNumContacts().
  ///
  std::vector<Vector6d> mu;

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
  /// SplitSolution::dimf().
  /// @return Reference to the stack of the Lagrange multipliers w.r.t. the
  /// switching constraints.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> xi_stack();

  ///
  /// @brief const version of SplitSolution::xi_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> xi_stack() const;

  ///
  /// @brief Returns the dimension of the contact at the current contact status.
  /// @return Dimension of the contact.
  ///
  int dimf() const;

  ///
  /// @brief Returns the dimension of the switching constraint.
  /// @return Dimension of the switching constraint.
  ///
  int dims() const;

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
  /// @param[in] impulse Flaf if this is impulse stage or not.
  ///
  void integrate(const Robot& robot, const double step_size, 
                 const SplitDirection& d, const bool impulse=false);

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
  /// @brief Return L-infinity Norm of the lagrange multipliers. Used in
  /// line search.
  ///
  double lagrangeMultiplierLinfNorm() const;

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
  std::vector<ContactType> contact_types_;
  std::vector<bool> is_contact_active_;
  int dimf_, dims_, max_num_contacts_;

};

} // namespace robotoc 

#include "robotoc/core/split_solution.hxx"

#endif // ROBOTOC_SPLIT_SOLUTION_HPP_