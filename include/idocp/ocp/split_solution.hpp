#ifndef IDOCP_SPLIT_SOLUTION_HPP_
#define IDOCP_SPLIT_SOLUTION_HPP_

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
      mu(Eigen::VectorXd::Zero(robot.dim_passive()+robot.max_dimf())),
      a(Eigen::VectorXd::Zero(robot.dimv())),
      f(Eigen::VectorXd::Zero(robot.max_dimf())),
      q(Eigen::VectorXd::Zero(robot.dimq())),
      v(Eigen::VectorXd::Zero(robot.dimv())),
      u(Eigen::VectorXd::Zero(robot.dimv())),
      beta(Eigen::VectorXd::Zero(robot.dimv())),
      dimc_(robot.dim_passive()+robot.dimf()),
      dimf_(robot.dimf()) {
    robot.normalizeConfiguration(q);
  }

  SplitSolution() 
    : lmd(),
      gmm(),
      mu(),
      a(),
      f(),
      q(),
      v(),
      u(),
      beta(),
      dimc_(0),
      dimf_(0) {
  }

  ~SplitSolution() {
  }

  SplitSolution(const SplitSolution&) = default;

  SplitSolution& operator=(const SplitSolution&) = default;
 
  SplitSolution(SplitSolution&&) noexcept = default;

  SplitSolution& operator=(SplitSolution&&) noexcept = default;

  inline void setContactStatus(const Robot& robot) {
    dimc_ = robot.dim_passive() + robot.dimf();
    dimf_ = robot.dimf();
  }

  inline Eigen::Ref<Eigen::VectorXd> f_active() {
    return f.head(dimf_);
  }

  inline Eigen::Ref<Eigen::VectorXd> mu_active() {
    return mu.head(dimc_);
  }

  inline Eigen::Ref<const Eigen::VectorXd> f_active() const {
    return f.head(dimf_);
  }

  inline Eigen::Ref<const Eigen::VectorXd> mu_active() const {
    return mu.head(dimc_);
  }

  int dimc() const {
    return dimc_;
  }

  int dimf() const {
    return dimf_;
  }

  Eigen::VectorXd lmd, gmm, mu, a, f, q, v, u, beta;

private:
  int dimc_, dimf_;

};

} // namespace idocp 


#endif // IDOCP_SPLIT_SOLUTION_HPP_