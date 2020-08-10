#ifndef IDOCP_SPLIT_OCP_SOLUTION_HPP_
#define IDOCP_SPLIT_OCP_SOLUTION_HPP_

#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class SplitOCPSolution {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  SplitOCPSolution(const Robot& robot) 
    : lmd(Eigen::VectorXd::Zero(robot.dimv())),
      gmm(Eigen::VectorXd::Zero(robot.dimv())),
      mu(Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf())),
      a(Eigen::VectorXd::Zero(robot.dimv())),
      f(Eigen::VectorXd::Zero(robot.max_dimf())),
      q(Eigen::VectorXd::Zero(robot.dimq())),
      v(Eigen::VectorXd::Zero(robot.dimv())),
      u(Eigen::VectorXd::Zero(robot.dimv())),
      beta(Eigen::VectorXd::Zero(robot.dimv())) {
  }

  SplitOCPSolution() 
    : lmd(),
      gmm(),
      mu(),
      a(),
      f(),
      q(),
      v(),
      u(),
      beta() {
  }

  ~SplitOCPSolution() {
  }

  SplitOCPSolution(const SplitOCPSolution&) = default;

  SplitOCPSolution& operator=(const SplitOCPSolution&) = default;
 
  SplitOCPSolution(SplitOCPSolution&&) noexcept = default;

  SplitOCPSolution& operator=(SplitOCPSolution&&) noexcept = default;

  Eigen::VectorXd lmd, gmm, mu, a, f, q, v, u, beta; 

private:
};

} // namespace idocp 


#endif // IDOCP_SPLIT_OCP_SOLUTION_HPP_