#ifndef INVDYNOCP_OCP_HPP_
#define INVDYNOCP_OCP_HPP_ 

#include <vector>
#include <omp.h>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/split_ocp.hpp"
#include "ocp/cost_function_interface.hpp"
#include "ocp/constraints_interface.hpp"


namespace invdynocp {

class OCP {
public:
  // Constructor. 
  OCP(const Robot& robot, const CostFunctionInterface* cost,
      const ConstraintsInterface* constraints, const double T, 
      const unsigned int N, const unsigned int num_proc);

  void solveSQP(const double t, const Eigen::VectorXd& x);

  // Prohibits copy constructor.
  OCP(const OCP&) = delete;

  // Prohibits copy operator.
  OCP& operator=(const OCP&) = delete;

private:
  std::vector<SplitOCP> split_ocps_;
  std::vector<Robot> robots_;
  CostFunctionInterface* cost_;
  ConstraintsInterface* constraints_;
  double T_, dtau_;
  unsigned int N_, num_proc_;
  std::vector<Eigen::VectorXd> x_, a_, lmd_, dx_, da_, dlmd_, s_;
  std::vector<Eigen::MatrixXd> P_;

};

} // namespace invdynocp


#endif // INVDYNOCP_OCP_HPP_