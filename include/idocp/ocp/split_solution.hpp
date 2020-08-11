#ifndef IDOCP_SPLIT_SOLUTION_HPP_
#define IDOCP_SPLIT_SOLUTION_HPP_


#include <assert.h>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/ocp/kkt_composition.hpp"


namespace idocp {

class SplitSolution {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  SplitSolution(const Robot& robot) 
    : lmd(Eigen::VectorXd::Zero(robot.dimv())),
      gmm(Eigen::VectorXd::Zero(robot.dimv())),
      a(Eigen::VectorXd::Zero(robot.dimv())),
      q(Eigen::VectorXd::Zero(robot.dimq())),
      v(Eigen::VectorXd::Zero(robot.dimv())),
      u(Eigen::VectorXd::Zero(robot.dimv())),
      beta(Eigen::VectorXd::Zero(robot.dimv())),
      mu_(Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf())),
      f_(Eigen::VectorXd::Zero(robot.max_dimf())),
      mu(mu_.head(robot.dim_passive()+robot.dimf())),
      f(f_.head(robot.dimf())) {
    robot.normalizeConfiguration(q);
  }

  SplitSolution() 
    : lmd(),
      gmm(),
      a(),
      q(),
      v(),
      u(),
      beta(),
      mu_(),
      f_(),
      mu(mu_),
      f(f_) {
  }

  ~SplitSolution() {
  }

  SplitSolution(const SplitSolution&) = default;

  SplitSolution& operator=(const SplitSolution&) = default;
 
  SplitSolution(SplitSolution&&) noexcept = default;

  SplitSolution& operator=(SplitSolution&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    mu = mu_.head(robot.dim_passive()+robot.dimf());
    f = f_.head(robot.dimf());
  }

  Eigen::VectorXd lmd, gmm, a, q, v, u, beta;
  Eigen::Ref<Eigen::VectorXd> mu, f;

private:
  Eigen::VectorXd mu_, f_;

};

} // namespace idocp 


#endif // IDOCP_SPLIT_SOLUTION_HPP_