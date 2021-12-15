#ifndef ROBOTOC_SPLIT_DIRECTION_HPP_
#define ROBOTOC_SPLIT_DIRECTION_HPP_

#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/contact_status.hpp"
#include "robotoc/robot/impulse_status.hpp"


namespace robotoc {

///
/// @class SplitDirection
/// @brief Newton direction of the solution to the optimal control problem 
/// split into a time stage. 
///
class SplitDirection {
public:
  ///
  /// @brief Construct a split solution.
  /// @param[in] robot Robot model. 
  ///
  SplitDirection(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  SplitDirection();

  ///
  /// @brief Destructor. 
  ///
  ~SplitDirection();

  ///
  /// @brief Default copy constructor. 
  ///
  SplitDirection(const SplitDirection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  SplitDirection& operator=(const SplitDirection&) = default;
 
  ///
  /// @brief Default move constructor. 
  ///
  SplitDirection(SplitDirection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  SplitDirection& operator=(SplitDirection&&) noexcept = default;

  ///
  /// @brief Sets contact status, i.e., set dimension of the contact forces.
  /// @param[in] contact_status Contact status.
  ///
  void setContactStatus(const ContactStatus& contact_status);

  ///
  /// @brief Sets the impulse status, i.e., set the dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Stack of the Newton directions of SplitSolution::q and 
  /// SplitSolution::v. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd dx;

  ///
  /// @brief Newton direction of SplitSolution::q. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dq();

  ///
  /// @brief const version of SplitDirection::dq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dq() const;

  ///
  /// @brief Newton direction of SplitSolution::gmm. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dv();

  ///
  /// @brief const version of SplitDirection::dv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dv() const;

  ///
  /// @brief Newton direction of SplitSolution::u. Size is Robot::dimu().
  ///
  Eigen::VectorXd du;

  ///
  /// @brief Stack of Newton direction of SplitSolution::a and SplitSolution::f. 
  /// Size is Robot::dimv() + ContactStatus::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> daf();

  ///
  /// @brief const version of SplitDirection::daf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> daf() const;

  ///
  /// @brief Newton direction of SplitSolution::a. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> da();

  ///
  /// @brief const version of SplitDirection::da().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> da() const;

  ///
  /// @brief Newton direction of SplitSolution::f_stack(). Size is 
  /// ContactStatus::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> df();

  ///
  /// @brief const version of SplitDirection::df().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> df() const;

  ///
  /// @brief Stack of the Newton direction of SplitSolution::lmd and 
  /// SplitSolution::gmm. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd dlmdgmm;

  ///
  /// @brief Newton direction of SplitSolution::lmd. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dlmd();

  ///
  /// @brief const version of SplitDirection::dlmd().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dlmd() const;

  ///
  /// @brief Newton direction of SplitSolution::gmm. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dgmm();

  ///
  /// @brief const version of SplitDirection::dgmm().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dgmm() const;

  ///
  /// @brief Stack of the Newton direction of SplitSolution::beta and 
  /// SplitSolution::mu_stack(). Size is Robot::dimv() + SplitSolution::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbetamu();

  ///
  /// @brief const version of SplitDirection::dbetamu(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbetamu() const;

  ///
  /// @brief Newton direction of SplitSolution::beta. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbeta();

  ///
  /// @brief const version of SplitDirection::dbeta().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbeta() const;

  ///
  /// @brief Newton direction of SplitSolution::mu_stack(). Size is 
  /// SplitSolution::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dmu();

  ///
  /// @brief const version of SplitDirection::dmu().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dmu() const;

  ///
  /// @brief Newton direction of SplitSolution::nu_passive. Size is 
  /// Robot::dim_passive().
  ///
  Eigen::VectorXd dnu_passive;

  ///
  /// @brief Newton direction of SplitSolution::xi_stack(). Size is 
  /// SplitSolution::dimi().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dxi();

  ///
  /// @brief const version of SplitDirection::dxi().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dxi() const;

  ///
  /// @brief Newton direction of the switching time.
  ///
  double dts;

  ///
  /// @brief Newton direction of the next switching time.
  ///
  double dts_next;

  ///
  /// @brief Set the all directions zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the contact.
  /// @return Dimension of the contact.
  ///
  int dimf() const;

  ///
  /// @brief Returns the dimension of the impulse.
  /// @return Dimension of the impulses.
  ///
  int dimi() const;

  ///
  /// @brief Checks dimensional consistency of each component. 
  /// @return true if the dimension is consistent. false if not.
  ///
  bool isDimensionConsistent() const;

  ///
  /// @brief Return true if two SplitDirection have the same values and false if 
  /// not. 
  /// @param[in] other Split direction that is compared with this object.
  ///
  bool isApprox(const SplitDirection& other) const;

  ///
  /// @brief Sets each component vector by random value based on the current 
  /// contact status. 
  ///
  void setRandom();

  ///
  /// @brief Sets each component vector by random value. Contact status is reset.
  /// @param[in] contact_status Contact status.
  ///
  void setRandom(const ContactStatus& contact_status);

  ///
  /// @brief Sets each component vector by random value. Impulse status is reset.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const ImpulseStatus& impulse_status);

  ///
  /// @brief Sets each component vector by random value. Contact status and 
  /// impulse status are reset.
  /// @param[in] contact_status Contact status.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const ContactStatus& contact_status,
                 const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  static SplitDirection Random(const Robot& robot);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status.
  ///
  static SplitDirection Random(const Robot& robot, 
                               const ContactStatus& contact_status);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] impulse_status Impulse status.
  ///
  static SplitDirection Random(const Robot& robot, 
                               const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] contact_status Contact status.
  /// @param[in] impulse_status Impulse status.
  ///
  static SplitDirection Random(const Robot& robot, 
                               const ContactStatus& contact_status,
                               const ImpulseStatus& impulse_status);

  ///
  /// @brief Displays the split direction onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, const SplitDirection& d);

private:
  Eigen::VectorXd daf_full_, dbetamu_full_, dxi_full_;
  int dimv_, dimu_, dim_passive_, dimf_, dimi_;

};

} // namespace robotoc 

#include "robotoc/ocp/split_direction.hxx"

#endif // ROBOTOC_SPLIT_OCP_DIRECTION_HPP_