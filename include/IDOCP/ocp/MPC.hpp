#ifndef IDOCP_MPC_HPP_
#define IDOCP_MPC_HPP_ 

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/ocp.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"


namespace idocp {

class MPC {
public:
  // Constructor. 
  MPC(const Robot& robot, const CostFunctionInterface* cost,
      const ConstraintsInterface* constraints, const double T, 
      const unsigned int N, const unsigned int num_proc);

  void initializeSolution(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, 
                          const unsigned int max_itr=100);

  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  void getControlInput(Eigen::VectorXd& u);

  double optimalityError(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v);

  // Prohibits copy constructor.
  MPC(const MPC&) = delete;

  // Prohibits copy operator.
  MPC& operator=(const MPC&) = delete;

private:
  OCP ocp_;
};

} // namespace idocp 


#endif // IDOCP_MPC_HPP_