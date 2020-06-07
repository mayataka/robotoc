#ifndef INVDYNOCP_SPLIT_OCP_HPP_
#define INVDYNOCP_SPLIT_OCP_HPP_

#include <vector>

#include "Eigen/Core"

#include "robot/robot.hpp"
#include "ocp/cost_function_interface.hpp"
#include "ocp/constraints_interface.hpp"


namespace invdynocp {

class SplitOCP {
public:
  // Constructor. Does not allocate raw arrays.
  // Argments:
  //    model: The pinocchio model. Before call this function, pinocchio model
  //      must be initialized by pinocchio::buildModel().
  SplitOCP(const Robot& robot, const CostFunctionInterface* cost,
           const ConstraintsInterface* constraints);

  ~SplitOCP();
 
  // Copy constructor.
  SplitOCP(const SplitOCP& other) = default;

  // Copy operator.
  SplitOCP& operator=(const SplitOCP& other) = default;

  void linearizeOCP(Robot& robot, const double t, const double dtau, 
                    const Eigen::VectorXd& lmd, const Eigen::VectorXd& x, 
                    const Eigen::VectorXd& a, const Eigen::VectorXd& lmd_next, 
                    const Eigen::VectorXd& x_next);

  void linearizeTerminalCost(Robot& robot, const double t, 
                             const Eigen::VectorXd& x, Eigen::MatrixXd& Qxx, 
                             Eigen::VectorXd& Qx);

  void backwardRecursion(const double dtau, const Eigen::MatrixXd& P1, 
                         const Eigen::VectorXd& s1, Eigen::MatrixXd& P, 
                         Eigen::VectorXd& s);

  void forwardRecursion(const double dtau, const Eigen::VectorXd& dx, 
                        Eigen::VectorXd& da, Eigen::VectorXd& dx_next) const;

  void updateOCP(Robot& robot, const Eigen::VectorXd& dx, 
                 const Eigen::VectorXd& da, const Eigen::MatrixXd& P, 
                 const Eigen::VectorXd& s, Eigen::VectorXd& x, 
                 Eigen::VectorXd& a, Eigen::VectorXd& lmd);

private:
  CostFunctionInterface *cost_;
  ConstraintsInterface *constraints_;
  unsigned int dimq_, dimv_;
  Eigen::VectorXd u_, lu_, lx_, la_, k_;
  Eigen::VectorXd x_res_, lmd_res_, a_res_;
  Eigen::MatrixXd du_dx_, du_da_, Qxx_, Qxa_, Qaa_, Ginv_, K_;

};

} // namespace invdynocp


#endif // INVDYNOCP_SPLIT_OCP_HPP_