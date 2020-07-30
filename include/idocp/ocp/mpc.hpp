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
      const int N, const int num_proc=1);

  ~MPC();

  // Use default copy constructor.
  MPC(const MPC&) = default;

  // Use default copy operator.
  MPC& operator=(const MPC&) = default;

  void initializeSolution(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, 
                          const int max_itr=100);

  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v);

  void getControlInput(Eigen::VectorXd& u);

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv);

  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v);

private:
  OCP ocp_;
};

} // namespace idocp 


#endif // IDOCP_MPC_HPP_