#ifndef IDOCP_RICCATI_CONSTRAINT_SOLUTION_HPP_
#define IDOCP_RICCATI_CONSTRAINT_SOLUTION_HPP_

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

class RiccatiConstraintSolution {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  RiccatiConstraintSolution(const Robot& robot, const int N, 
                            const int max_impulse_stage);

  // Default constructor.
  RiccatiConstraintSolution();

  // Destructor.
  ~RiccatiConstraintSolution();
 
  // Default copy constructor.
  RiccatiConstraintSolution(const RiccatiConstraintSolution&) = default;

  // Default copy operator.
  RiccatiConstraintSolution& operator=(const RiccatiConstraintSolution&) 
      = default;

  // Default move constructor.
  RiccatiConstraintSolution(RiccatiConstraintSolution&&) noexcept = default;

  // Default move assign operator.
  RiccatiConstraintSolution& operator=(RiccatiConstraintSolution&&) noexcept 
      = default;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

  bool isActive() const;

  Eigen::Block<Eigen::MatrixXd> T(const int i);

  const Eigen::Block<const Eigen::MatrixXd> T(const int i) const;

  Eigen::Block<Eigen::MatrixXd> T_impulse(const int i);

  const Eigen::Block<const Eigen::MatrixXd> T_impulse(const int i) const;

  Eigen::Block<Eigen::MatrixXd> ENEt();

  const Eigen::Block<const Eigen::MatrixXd> ENEt() const;

  Eigen::Block<Eigen::MatrixXd> EqNqq();

  const Eigen::Block<const Eigen::MatrixXd> EqNqq() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::MatrixXd> T_full_, T_impulse_full_;
  Eigen::MatrixXd EqNqq_full_, ENEt_full_;
  Eigen::LLT<Eigen::MatrixXd> llt_;
  int max_impulse_stage_, dimv_, dimx_, dimf_;
  bool is_active_;

};

} // namespace idocp

#include "idocp/ocp/riccati_constraint_solution.hxx"

#endif // IDOCP_RICCATI_CONSTRAINT_SOLUTION_HPP_ 