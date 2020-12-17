#ifndef IDOCP_SPLIT_DIRECTION_HPP_
#define IDOCP_SPLIT_DIRECTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/contact_status.hpp"


namespace idocp {

///
/// @class SplitDirection
/// @brief Newton direction of SplitSolution. 
///
class SplitDirection {
public:

  using Vector6d = Eigen::Matrix<double, 6, 1>;

  ///
  /// @brief Construct a split solution.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
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
  /// @brief Set contact status, i.e., set dimension of the contact.
  /// @param[in] dimf Total dimension of the contact.
  ///
  void setContactStatusByDimension(const int dimf);

  ///
  /// @brief Returns the Newton direction of SplitSolution::lmd and 
  /// SplitSolution::gmm. Size is 2 * Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::lmd and 
  /// SplitSolution::gmm.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dlmdgmm();

  ///
  /// @brief Const verison of SplitDirection::dlmdgmm().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dlmdgmm() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::lmd. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::lmd.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dlmd();

  ///
  /// @brief Const version of SplitDirection::dlmd().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dlmd() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::gmm. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::gmm.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dgmm();

  ///
  /// @brief Const version of SplitDirection::dgmm().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dgmm() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::u. Size is 
  /// Robot::dimu().
  /// @return Reference to the Newton direction of SplitSolution::u.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> du();

  ///
  /// @brief Const version of SplitDirection::du().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> du() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::q. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::q.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dq();

  ///
  /// @brief Const version of SplitDirection::dq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dq() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::gmm. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::gmm.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dv();

  ///
  /// @brief Const version of SplitDirection::dv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dv() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::q and 
  /// SplitSolution::v. Size is 2 * Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::q and 
  /// SplitSolution::v.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dx();

  ///
  /// @brief Const version of SplitDirection::dx().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dx() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::a and 
  /// SplitSolution::f. Size is ContactStatus::dimf() + Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::a and 
  /// SplitSolution::f.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> daf();

  ///
  /// @brief Const version of SplitDirection::daf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> daf() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::a. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::a.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> da();

  ///
  /// @brief Const version of SplitDirection::da().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> da() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::f_stack(). Size is 
  /// ContactStatus::dimf().
  /// @return Reference to the Newton direction of SplitSolution::f_stack().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> df();

  ///
  /// @brief Const version of SplitDirection::df().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> df() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::beta and 
  /// SplitSolution::mu_stack(). Size is Robot::dimv() + SplitSolution::dimf().
  /// @return Reference to the Newton direction of SplitSolution::beta and 
  /// SplitSolution::mu_stack().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbetamu();

  ///
  /// @brief Const version of SplitDirection::dbetamu(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbetamu() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::beta. Size 
  /// is Robot::dimv().
  /// @return Reference to the Newton direction of SplitSolution::beta.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbeta();

  ///
  /// @brief Const version of SplitDirection::dbeta().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbeta() const;

  ///
  /// @brief Returns the Newton direction of SplitSolution::mu_stack(). Size is 
  /// SplitSolution::dimf().
  /// @return Reference to the Newton direction of SplitSolution::mu_stack().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dmu();

  ///
  /// @brief Const version of SplitDirection::dmu().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dmu() const;

  ///
  /// @brief Set the all directions zero.
  ///
  void setZero();

  ///
  /// @brief Returns the dimension of the condensed KKT conditions.
  /// @return Dimension of the condensed KKT conditions.
  ///
  int dimKKT() const;

  ///
  /// @brief Returns the dimension of the stack of contact forces at the current 
  /// contact status.
  /// @return Dimension of the stack of contact forces.
  ///
  int dimf() const;

  /// @brief Stack of the Newton direction composed of SplitDirection::dlmd(), 
  /// SplitDirection::dgmm(), SplitDirection::du(), SplitDirection::dq(), 
  /// and SplitDirection::dv(). Size is 4 * Robot::dimv() + Robot::dimu().
  Eigen::VectorXd split_direction;

  /// @brief Newton direction of SplitSolution::u_passive. Size is 6.
  Vector6d du_passive;

  /// @brief Newton direction of SplitSolution::nu_passive. Size is 6.
  Vector6d dnu_passive;

  ///
  /// @brief Return true if two SplitDirection have the same values and false if 
  /// not. 
  /// @param[in] other Split direction that is compared with this object.
  ///
  bool isApprox(const SplitDirection& other) const;

  ///
  /// @brief Set each component vector by random value based on the current 
  /// contact status. 
  ///
  void setRandom();

  ///
  /// @brief Set each component vector by random value. Contact status is reset.
  /// @param[in] contact_status Contact status.
  ///
  void setRandom(const ContactStatus& contact_status);

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

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd daf_full_, dbetamu_full_;
  int dimv_, dimu_, dimx_, dimf_, dimKKT_;
  bool has_floating_base_;

};

} // namespace idocp 

#include "idocp/ocp/split_direction.hxx"

#endif // IDOCP_SPLIT_OCP_DIRECTION_HPP_