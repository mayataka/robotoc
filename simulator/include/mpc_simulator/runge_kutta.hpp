#ifndef IDOCP_SIMULATOR_RUNGE_KUTTA_HPP_
#define IDOCP_SIMULATOR_RUNGE_KUTTA_HPP_

#include "Eigen/Core"
#include "robot/robot.hpp"


namespace idocp {
namespace simulator {

class RungeKutta {
public:
  RungeKutta(const Robot& robot);

  // Prohibit copy constructor.
  RungeKutta(const RungeKutta& other) = delete;

  // Prohibit copy operator.
  RungeKutta& operator=(const RungeKutta& other) = delete;

  void integrate(const double integration_length, const Eigen::VectorXd& q, 
                 const Eigen::VectorXd& v, const Eigen::VectorXd& tau,
                 Eigen::VectorXd& q_next, Eigen::VectorXd& v_next);

private:
  Robot robot_;
  Eigen::VectorXd dq_, dv_, kq1_, kq2_, kq3_, kq4_, kv1_, kv2_, kv3_, kv4_;

};

} // namespace simulator
} // namespace idocp


#endif // IDOCP_SIMULATOR_RUNGE_KUTTA_HPP_ 