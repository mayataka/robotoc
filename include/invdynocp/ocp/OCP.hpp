#ifndef INVDYNOCP_OCP_HPP_
#define INVDYNOCP_OCP_HPP_ 

#include <vector>
#include <omp.h>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/split_ocp.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"


namespace invdynocp {

class OCP {
public:
  // Constructor. 
  OCP(const Robot& robot, const CostFunctionInterface* cost,
      const ConstraintsInterface* constraints, const double T, 
      const unsigned int N, const unsigned int num_proc);

  void solveSQP(const double t, const Eigen::VectorXd& q, 
                const Eigen::VectorXd& v);

  void getInitialControlInput(Eigen::VectorXd& u);

  void setStateTrajectory(Eigen::VectorXd& q0, Eigen::VectorXd& v0, 
                          Eigen::VectorXd& qN, Eigen::VectorXd& vN);

  void printSolution() const;

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

};

} // namespace invdynocp


#endif // INVDYNOCP_OCP_HPP_