#ifndef IDOCP_SPLIT_SOLUTION_HPP_
#define IDOCP_SPLIT_SOLUTION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"
#include "idocp/ocp/split_direction.hpp"


namespace idocp {

///
/// @class SplitSolution
/// @brief Solution of the optimal control problem split into each time stage. 
///
class SplitSolution {
public:

  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Construct a split solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
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
  /// @brief Set contact status, i.e., set dimension of the contacts.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Stack of active contact forces. Size is ContactStatus::dimf().
  /// @return Reference to the stack of active contact forces.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> f_stack();

  ///
  /// @brief Const version of SplitSolution::f_stack().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> f_stack() const;

  ///
  /// @brief Set SplitSolution::f_stack() from SplitSolution::f.
  ///
  void set_f_stack();

  ///
  /// @brief Set SplitSolution::f from SplitSolution::f_stack().
  ///
  void set_f_vector();

  ///
  /// @brief Stack of Lagrange multiplier with respect to contact constraints 
  /// that is active at the current contact status. Size is 
  /// SplitSolution::dimf().
  /// @return Reference to the stack of Lagrange multiplier with respect to 
  /// contact constraints.
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
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of contact forces.
  ///
  int dimf() const;

  ///
  /// @brief Lagrange multiplier with respect to transition of SplitSolution::q. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd lmd;

  ///
  /// @brief Lagrange multiplier with respect to transition of SplitSolution::v. 
  /// Size is Robot::dimv().
  ///
  Eigen::VectorXd gmm;

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
  /// @brief Contact forces. 
  /// Size is Robot::maxPointContacts().
  ///
  std::vector<Eigen::Vector3d> f;

  ///
  /// @brief Control input torques. Size is Robot::dimu().
  ///
  Eigen::VectorXd u;

  ///
  /// @brief Lagrange multiplier with respect to inverse dynamics. Size is 
  /// Robot::dimv().
  ///
  Eigen::VectorXd beta;

  ///
  /// @brief Lagrange multiplier with respect to contact constraint. 
  /// Size is Robot::maxPointContacts().
  ///
  std::vector<Eigen::Vector3d> mu;

  ///
  /// @brief Control input torques of the virtual floating base joint. Size is 6.
  ///
  Vector6d u_passive;

  ///
  /// @brief Lagrange multiplier with respect to floating base. Size is 6.
  ///
  Vector6d nu_passive;

  ///
  /// @brief Return true if a contact is active and false if not.
  /// @param[in] contact_index Index of a contact of interedted. 
  /// @return true if a contact is active and false if not. 
  ///
  bool isContactActive(const int contact_index) const;

  ///
  /// @brief Return true if there are active contacts and false if not.
  /// @return true if there are active contacts and false if not. 
  ///
  bool hasActiveContacts() const;

  ///
  /// @brief Integrates the solution based on step size and direction. 
  /// @param[in] robot Robot model.
  /// @param[in] step_size Step size.
  /// @param[in] d Split direction.
  ///
  void integrate(const Robot& robot, const double step_size, 
                 const SplitDirection& d);

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

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd mu_stack_, f_stack_;
  bool has_floating_base_, has_active_contacts_;
  std::vector<bool> is_contact_active_;
  int dimf_;

};

} // namespace idocp 

#include "idocp/ocp/split_solution.hxx"

#endif // IDOCP_SPLIT_SOLUTION_HPP_