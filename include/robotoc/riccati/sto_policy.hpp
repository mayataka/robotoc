#ifndef ROBOTOC_STO_POLICY_HPP_
#define ROBOTOC_STO_POLICY_HPP_

#include "Eigen/Core"

#include "robotoc/robot/robot.hpp"


namespace robotoc {

///
/// @class STOPolicy
/// @brief The state feedback and feedforward policy of the switching time 
/// optimization (STO).
///
class STOPolicy {
public:
  ///
  /// @brief Constructs STO gain and feedforward term.
  /// @param[in] robot Robot model. 
  ///
  STOPolicy(const Robot& robot)
    : dtsdx(Eigen::VectorXd::Zero(2*robot.dimv())),
      dtsdts(0),
      dts0(0) {
  }

  ///
  /// @brief Default constructor. 
  ///
  STOPolicy() 
    : dtsdx(),
      dtsdts(0),
      dts0(0) {
  }

  ///
  /// @brief Destructor. 
  ///
  ~STOPolicy() {
  }

  ///
  /// @brief Default copy constructor. 
  ///
  STOPolicy(const STOPolicy&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  STOPolicy& operator=(const STOPolicy&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  STOPolicy(STOPolicy&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  STOPolicy& operator=(STOPolicy&&) noexcept = default;

  ///
  /// @brief Feedback gain of the STO policy, i.e. the sensitivity w.r.t. the 
  /// state. Size is 2 * Robot::dimv().
  ///
  Eigen::VectorXd dtsdx;

  ///
  /// @brief Feedback gain of the STO policy, i.e. the sensitivity w.r.t. the 
  /// previous switching time.
  ///
  double dtsdts;

  ///
  /// @brief Feedforward term of the STO policy.
  ///
  double dts0;

private:

};

} // namespace robotoc 

#endif // ROBOTOC_STO_POLICY_HPP_ 