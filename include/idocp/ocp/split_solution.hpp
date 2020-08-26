#ifndef IDOCP_SPLIT_SOLUTION_HPP_
#define IDOCP_SPLIT_SOLUTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"


namespace idocp {

class SplitSolution {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  SplitSolution(const Robot& robot);

  SplitSolution();

  ~SplitSolution();

  SplitSolution(const SplitSolution&) = default;

  SplitSolution& operator=(const SplitSolution&) = default;
 
  SplitSolution(SplitSolution&&) noexcept = default;

  SplitSolution& operator=(SplitSolution&&) noexcept = default;

  void setContactStatus(const Robot& robot);

  Eigen::VectorBlock<Eigen::VectorXd> f_active();

  Eigen::VectorBlock<Eigen::VectorXd> mu_active();

  const Eigen::VectorBlock<const Eigen::VectorXd> f_active() const;

  const Eigen::VectorBlock<const Eigen::VectorXd> mu_active() const;

  int dimc() const;

  int dimf() const;

  Eigen::VectorXd lmd, gmm, mu, a, f, q, v, u, beta;

private:
  int dimc_, dimf_;

};

} // namespace idocp 

#include "idocp/ocp/split_solution.hxx"

#endif // IDOCP_SPLIT_SOLUTION_HPP_