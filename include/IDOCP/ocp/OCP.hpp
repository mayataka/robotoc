#ifndef IDOCP_OCP_HPP_
#define IDOCP_OCP_HPP_ 

#include <vector>
#include <utility>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/split_OCP.hpp"
#include "ocp/line_search_filter.hpp"
#include "cost/cost_function_interface.hpp"
#include "constraints/constraints_interface.hpp"


namespace idocp {

class OCP {
public:
  // Constructor. 
  OCP(const Robot& robot, const CostFunctionInterface* cost,
      const ConstraintsInterface* constraints, const double T, 
      const unsigned int N, const unsigned int num_proc);

  void solveSQP(const double t, const Eigen::VectorXd& q, 
                const Eigen::VectorXd& v, bool use_line_search=true);

  void getInitialControlInput(Eigen::VectorXd& u);

  void setStateTrajectory(const Eigen::VectorXd& q, const Eigen::VectorXd& v);

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
  LineSearchFilter filter_;
  CostFunctionInterface* cost_;
  ConstraintsInterface* constraints_;
  double T_, dtau_, step_size_reduction_rate_, min_step_size_;
  unsigned int N_, num_proc_;
  std::vector<Eigen::VectorXd> q_, v_, a_, u_, beta_, lmd_, gmm_, 
                               dq_, dv_, da_, dlmd_, dgmm_, sq_, sv_;
  std::vector<Eigen::MatrixXd> Pqq_, Pqv_, Pvq_, Pvv_;
  Eigen::VectorXd primal_step_sizes_, dual_step_sizes_, cost_origin_, 
                  cost_search_, constraints_residual_origin_, 
                  constraints_residual_search_;

};

} // namespace idocp 


#endif // IDOCP_OCP_HPP_