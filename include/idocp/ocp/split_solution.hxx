#ifndef IDOCP_SPLIT_SOLUTION_HXX_
#define IDOCP_SPLIT_SOLUTION_HXX_

namespace idocp {

inline SplitSolution::SplitSolution(const Robot& robot) 
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


inline SplitSolution::SplitSolution() 
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


inline SplitSolution::~SplitSolution() {
}


inline void SplitSolution::setContactStatus(const Robot& robot) {
  dimc_ = robot.dim_passive() + robot.dimf();
  dimf_ = robot.dimf();
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::f_active() {
  return f.head(dimf_);
}


inline Eigen::VectorBlock<Eigen::VectorXd> SplitSolution::mu_active() {
  return mu.head(dimc_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::f_active() const {
  return f.head(dimf_);
}


inline const Eigen::VectorBlock<const Eigen::VectorXd> 
SplitSolution::mu_active() const {
  return mu.head(dimc_);
}


inline int SplitSolution::dimc() const {
  return dimc_;
}


inline int SplitSolution::dimf() const {
  return dimf_;
}

} // namespace idocp 

#endif // IDOCP_SPLIT_SOLUTION_HXX_