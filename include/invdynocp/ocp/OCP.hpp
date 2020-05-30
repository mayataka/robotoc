#ifndef INVDYNOCP_OCP_HPP_
#define INVDYNOCP_OCP_HPP_ 

#include <omp.h>
#include <string>
#include <vector>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/cost_function_interface.hpp"
#include "ocp/constraints_interface.hpp"


namespace invdynocp {

class OCP {
public:
  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized by pinocchio::buildModel().
  OCP(const Robot& robot, const CostFunctionInterface* cost,
      const ConstraintsInterface* constraints, const double T, 
      const unsigned int N, const unsigned int num_proc);

  // Destructor. 
  ~OCP();

  void solveOCP(const double t, const Eigen::VectorXd& x);

  void getInitialContorlInpnut(Eigen::VectorXd& u);

  // Returns the dimension of the torques correspoinding to the passive joints.
  unsigned int dimq() const;

  unsigned int dimv() const;

  unsigned int dimu() const;

  // Prohibits copy constructor.
  OCP(const OCP&) = default;

  // Prohibits copy operator.
  OCP& operator=(const OCP&) = default;

private:
  std::vector<Robot> robots_;
  CostFunctionInterface* cost_;
  ConstraintsInterface* constraints_;
  double T_;
  unsigned int N_, num_proc_, dimx_, dimTx_, dimu_, dimc_;

};

} // namespace invdynocp


#endif // INVDYNOCP_OCP_HPP_