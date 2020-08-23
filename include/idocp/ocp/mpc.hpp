#ifndef IDOCP_MPC_HPP_
#define IDOCP_MPC_HPP_ 

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/ocp.hpp"
#include "idocp/ocp/parnmpc.hpp"
#include "idocp/cost/cost_function.hpp"
#include "idocp/constraints/constraints.hpp"


namespace idocp {

template <typename OCPType>
class MPC {
public:
  // Constructor. 
  MPC(const Robot& robot, const std::shared_ptr<CostFunction>& cost,
      const std::shared_ptr<Constraints>& constraints, const double T, 
      const int N, const int num_proc=1)
    : ocp(robot, cost, constraints, T, N, num_proc) {
  }

  ~MPC() {}

  // Use default copy constructor.
  MPC(const MPC&) = default;

  // Use default copy operator.
  MPC& operator=(const MPC&) = default;

  // Use default move constructor.
  MPC(MPC&&) noexcept = default;

  // Use default move operator.
  MPC& operator=(MPC&&) noexcept = default;

  void initializeSolution(const double t, const Eigen::VectorXd& q, 
                          const Eigen::VectorXd& v, 
                          const int max_itr=100) {
    ocp_.setStateTrajectory(q, v,);
    for (int i=0; i<max_itr; ++i) {
      ocp_.updateSolution(t, q, v, true);
    }
  }

  void updateSolution(const double t, const Eigen::VectorXd& q, 
                      const Eigen::VectorXd& v) {
    ocp_.updateSolution(t, q, v, false);
  }

  void getControlInput(Eigen::VectorXd& u) {
    ocp_.getControlInput(u);
  }

  void getStateFeedbackGain(Eigen::MatrixXd& Kq, Eigen::MatrixXd& Kv);

  double KKTError(const double t, const Eigen::VectorXd& q, 
                  const Eigen::VectorXd& v) {
    return ocp_.KKTError(t, q, v);
  }

private:
  OCPType ocp;
};

} // namespace idocp 


#endif // IDOCP_MPC_HPP_