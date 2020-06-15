#ifndef IDOCP_OCP_HPP_
#define IDOCP_OCP_HPP_ 

#include <vector>
#include <omp.h>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/split_OCP.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"


namespace idocp {

class OCP {
public:
  // Constructor. 
  OCP(const Robot& robot, const CostFunctionInterface* cost,
      const ConstraintsInterface* constraints, const double T, 
      const unsigned int N, const unsigned int num_proc);

  void solveLQR(const double t, const Eigen::VectorXd& q, 
                const Eigen::VectorXd& v);

  void getInitialControlInput(Eigen::VectorXd& u);

  void setStateTrajectory(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

  void setStateTrajectory(const Eigen::VectorXd& q0, const Eigen::VectorXd& v0, 
                          const Eigen::VectorXd& qN, const Eigen::VectorXd& vN);

  double optimalityError(const double t, const Eigen::VectorXd& q, 
                         const Eigen::VectorXd& v);

  void printSolution();

  // Prohibits copy constructor.
  OCP(const OCP&) = delete;

  // Prohibits copy operator.
  OCP& operator=(const OCP&) = delete;

private:
  std::vector<SplitOCP> split_OCPs_;
  std::vector<Robot> robots_;
  CostFunctionInterface* cost_;
  ConstraintsInterface* constraints_;
  double T_, dtau_;
  unsigned int N_, num_proc_;
  std::vector<Eigen::VectorXd> q_, v_, a_, lmd_, gmm_, 
                               dq_, dv_, da_, dlmd_, dgmm_, sq_, sv_;
  std::vector<Eigen::MatrixXd> Pqq_, Pqv_, Pvq_, Pvv_;
  Eigen::VectorXd step_length_;

};

} // namespace idocp 


#endif // IDOCP_OCP_HPP_