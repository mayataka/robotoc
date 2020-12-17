#ifndef IDOCP_IMPULSE_SPLIT_DIRECTION_HPP_ 
#define IDOCP_IMPULSE_SPLIT_DIRECTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"


namespace idocp {

///
/// @class ImpulseSplitDirection
/// @brief Newton direction of an impulse stage. 
///
class ImpulseSplitDirection {
public:
  ///
  /// @brief Construct a impulse split direction.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  ImpulseSplitDirection(const Robot& robot);

  ///
  /// @brief Default constructor.  
  ///
  ImpulseSplitDirection();

  ///
  /// @brief Destructor. 
  ///
  ~ImpulseSplitDirection();

  ///
  /// @brief Use default copy constructor. 
  ///
  ImpulseSplitDirection(const ImpulseSplitDirection&) = default;

  ///
  /// @brief Use default copy assign operator. 
  ///
  ImpulseSplitDirection& operator=(const ImpulseSplitDirection&) = default;
 
  ///
  /// @brief Use default move constructor. 
  ///
  ImpulseSplitDirection(ImpulseSplitDirection&&) noexcept = default;

  ///
  /// @brief Use default move assign operator. 
  ///
  ImpulseSplitDirection& operator=(ImpulseSplitDirection&&) noexcept = default;

  ///
  /// @brief Set impulse status from robot model, i.e., set dimension of the 
  /// impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Set impulse status, i.e., set dimension of the impulse.
  /// @param[in] dimf Total dimension of the impulse.
  ///
  void setImpulseStatusByDimension(const int dimf);

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::lmd and 
  /// SplitSolution::gmm. Size is 2 * Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::lmd and 
  /// ImpulseSplitSolution::gmm.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dlmdgmm();

  ///
  /// @brief Const verison of ImpulseSplitDirection::dlmdgmm().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dlmdgmm() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::lmd. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::lmd.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dlmd();

  ///
  /// @brief const version of ImpulseSplitDirection::dlmd().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dlmd() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::gmm. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::gmm.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dgmm();

  ///
  /// @brief const version of ImpulseSplitDirection::dgmm().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dgmm() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::q. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::q.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dq();

  ///
  /// @brief const version of ImpulseSplitDirection::dq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dq() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::v. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::v.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dv();

  ///
  /// @brief const version of ImpulseSplitDirection::dv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dv() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::q and 
  /// ImpulseSplitSolution::v. Size is 2 * Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::q and 
  /// ImpulseSplitSolution::v.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dx();

  ///
  /// @brief const version of ImpulseSplitDirection::dx().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dx() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::dv and 
  /// ImpulseSplitSolution::f. Size is ImpulseStatus::dimf() + Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::a and 
  /// ImpulseSplitSolution::f.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> ddvf();

  ///
  /// @brief const version of ImpulseSplitDirection::ddvf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> ddvf() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::dv. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::dv.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> ddv();

  ///
  /// @brief Const version of ImpulseSplitDirection::ddv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> ddv() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::f_stack(). 
  /// Size is ImpulseStatus::dimf().
  /// @return Reference to the Newton direction of 
  /// ImpulseSplitSolution::f_stack().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> df();

  ///
  /// @brief Const version of ImpulseSplitDirection::df().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> df() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::beta and 
  /// ImpulseSplitSolution::mu_stack(). Size is Robot::dimv() + 
  /// ImpulseStatus::dimf().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::beta 
  /// and ImpulseSplitSolution::mu_stack().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbetamu();

  ///
  /// @brief const version of ImpulseSplitDirection::dbetamu(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbetamu() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::beta. Size 
  /// is Robot::dimv().
  /// @return Reference to the Newton direction of ImpulseSplitSolution::beta.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbeta();

  ///
  /// @brief Const version of ImpulseSplitDirection::dbeta().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbeta() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::mu_stack(). 
  /// Size is ImpulseSplitSolution::dimf().
  /// @return Reference to the Newton direction of 
  /// ImpulseSplitSolution::mu_stack().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dmu();

  ///
  /// @brief Const version of ImpulseSplitDirection::dmu().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dmu() const;

  ///
  /// @brief Returns the Newton direction of ImpulseSplitSolution::xi_stack(). 
  /// Size is ImpulseSplitSolution::dimf().
  /// @return Reference to the Newton direction of 
  /// ImpulseSplitSolution::xi_stack().
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dxi();

  ///
  /// @brief Const version of ImpulseSplitDirection::dxi().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dxi() const;

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
  /// SplitDirection::dgmm(), SplitDirection::dq(), 
  /// and SplitDirection::dv(). Size is 4 * Robot::dimv() + Robot::dimu().
  Eigen::VectorXd split_direction;

  ///
  /// @brief Return true if two SplitDirection have the same values and false if 
  /// not. 
  /// @param[in] other Split direction that is compared with this object.
  ///
  bool isApprox(const ImpulseSplitDirection& other) const;

  ///
  /// @brief Set each component vector by random value based on the current 
  /// contact status. 
  ///
  void setRandom();

  ///
  /// @brief Set each component vector by random value. Impulse status is reset.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  ///
  static ImpulseSplitDirection Random(const Robot& robot);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] impulse_status Impulse status.
  ///
  static ImpulseSplitDirection Random(const Robot& robot, 
                                      const ImpulseStatus& impulse_status);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  Eigen::VectorXd ddvf_full_, dbetamu_full_, dxi_full_;
  int dimv_, dimx_, dimf_, dimKKT_;
  bool has_floating_base_;

};

} // namespace idocp 

#include "idocp/impulse/impulse_split_direction.hxx"

#endif // IDOCP_IMPULSE_SPLIT_DIRECTION_HPP_ 