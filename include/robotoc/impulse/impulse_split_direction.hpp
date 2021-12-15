#ifndef ROBOTOC_IMPULSE_SPLIT_DIRECTION_HPP_ 
#define ROBOTOC_IMPULSE_SPLIT_DIRECTION_HPP_

#include <iostream>

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"
#include "robotoc/robot/impulse_status.hpp"


namespace robotoc {

///
/// @class ImpulseSplitDirection
/// @brief Newton direction of the solution to the optimal control problem 
/// split into an impulse stage. 
///
class ImpulseSplitDirection {
public:
  ///
  /// @brief Constructs a split direction.
  /// @param[in] robot Robot model. 
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
  /// @brief Default copy constructor. 
  ///
  ImpulseSplitDirection(const ImpulseSplitDirection&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  ImpulseSplitDirection& operator=(const ImpulseSplitDirection&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  ImpulseSplitDirection(ImpulseSplitDirection&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  ImpulseSplitDirection& operator=(ImpulseSplitDirection&&) noexcept = default;

  ///
  /// @brief Sets the impulse status, i.e., set the dimension of the impulse.
  /// @param[in] impulse_status Impulse status.
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief Stack of the Newton directions of ImpulseSplitSolution::q and 
  /// ImpulseSplitSolution::v. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd dx;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::q. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dq();

  ///
  /// @brief const version of ImpulseSplitDirection::dq().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dq() const;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::v. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dv();

  ///
  /// @brief const version of ImpulseSplitDirection::dv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dv() const;

  ///
  /// @brief Stack of the Newton directions of ImpulseSplitSolution::dv and 
  /// ImpulseSplitSolution::f. Size is Robot::dimv() + ImpulseStatus::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> ddvf();

  ///
  /// @brief const version of ImpulseSplitDirection::ddvf().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> ddvf() const;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::dv. Size is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> ddv();

  ///
  /// @brief const version of ImpulseSplitDirection::ddv().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> ddv() const;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::f_stack(). 
  /// Size is ImpulseStatus::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> df();

  ///
  /// @brief const version of ImpulseSplitDirection::df().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> df() const;

  ///
  /// @brief Stack of the Newton directions of ImpulseSplitSolution::lmd and 
  /// ImpulseSplitSolution::gmm. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd dlmdgmm;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::lmd. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dlmd();

  ///
  /// @brief const version of ImpulseSplitDirection::dlmd().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dlmd() const;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::gmm. Size is 
  /// Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dgmm();

  ///
  /// @brief const version of ImpulseSplitDirection::dgmm().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dgmm() const;

  ///
  /// @brief Stack of the Newton directions of ImpulseSplitSolution::beta and 
  /// ImpulseSplitSolution::mu_stack(). Size is Robot::dimv() + 
  /// ImpulseStatus::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbetamu();

  ///
  /// @brief const version of ImpulseSplitDirection::dbetamu(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbetamu() const;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::beta. Size 
  /// is Robot::dimv().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dbeta();

  ///
  /// @brief const version of ImpulseSplitDirection::dbeta().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dbeta() const;

  ///
  /// @brief Newton direction of ImpulseSplitSolution::mu_stack(). 
  /// Size is ImpulseSplitSolution::dimf().
  /// @return Reference to the Newton direction.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> dmu();

  ///
  /// @brief const version of ImpulseSplitDirection::dmu().
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> dmu() const;

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
  /// @brief Returns the dimension of the impulses.
  /// @return Dimension of the impulses.
  ///
  int dimi() const;

  ///
  /// @brief Checks dimensional consistency of each component. 
  /// @return true if the dimension is consistent. false if not.
  ///
  bool isDimensionConsistent() const;

  ///
  /// @brief Returns true if two ImpulseSplitDirection are the same and false
  /// if not. 
  /// @param[in] other Split impulse direction that is compared with this object.
  ///
  bool isApprox(const ImpulseSplitDirection& other) const;

  ///
  /// @brief Sets each component vector by random value based on the current 
  /// impulse status. 
  ///
  void setRandom();

  ///
  /// @brief Sets each component vector by random value. Impulse status is reset.
  /// @param[in] impulse_status Impulse status.
  ///
  void setRandom(const ImpulseStatus& impulse_status);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. 
  ///
  static ImpulseSplitDirection Random(const Robot& robot);

  ///
  /// @brief Generates split direction filled randomly.
  /// @return Split direction filled randomly.
  /// @param[in] robot Robot model. 
  /// @param[in] impulse_status Impulse status.
  ///
  static ImpulseSplitDirection Random(const Robot& robot, 
                                      const ImpulseStatus& impulse_status);

  ///
  /// @brief Displays the impulse split direction onto a ostream.
  ///
  void disp(std::ostream& os) const;

  friend std::ostream& operator<<(std::ostream& os, 
                                  const ImpulseSplitDirection& d);

private:
  Eigen::VectorXd ddvf_full_, dbetamu_full_;
  int dimv_, dimx_, dimi_;

};

} // namespace robotoc 

#include "robotoc/impulse/impulse_split_direction.hxx"

#endif // ROBOTOC_IMPULSE_SPLIT_DIRECTION_HPP_ 