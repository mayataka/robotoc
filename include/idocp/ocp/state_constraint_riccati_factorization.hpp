#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

class StateConstraintRiccatiFactorization {
public:
  // Constructor.
  // Argments:
  //    robot: The robot model that has been already initialized.
  StateConstraintRiccatiFactorization(const Robot& robot, const int N);

  // Default constructor.
  StateConstraintRiccatiFactorization();

  // Destructor.
  ~StateConstraintRiccatiFactorization();
 
  // Default copy constructor.
  StateConstraintRiccatiFactorization(
      const StateConstraintRiccatiFactorization&) = default;

  // Default copy operator.
  StateConstraintRiccatiFactorization& operator=(
      const StateConstraintRiccatiFactorization&) = default;

  // Default move constructor.
  StateConstraintRiccatiFactorization(
      StateConstraintRiccatiFactorization&&) noexcept = default;

  // Default move assign operator.
  StateConstraintRiccatiFactorization& operator=(
      StateConstraintRiccatiFactorization&&) noexcept = default;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

  Eigen::Block<Eigen::MatrixXd> T(const int time_stage);

  const Eigen::Block<const Eigen::MatrixXd> T(const int time_stage) const;

  Eigen::Block<Eigen::MatrixXd> T_aux(const int time_stage);

  const Eigen::Block<const Eigen::MatrixXd> T_aux(const int time_stage) const;

  Eigen::Block<Eigen::MatrixXd> ENEt();

  const Eigen::Block<const Eigen::MatrixXd> ENEt() const;

  Eigen::Block<Eigen::MatrixXd> EqNqq();

  const Eigen::Block<const Eigen::MatrixXd> EqNqq() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::MatrixXd> T_full_, T_aux_full_;
  Eigen::MatrixXd EqNqq_full_, ENEt_full_;
  int N_, dimv_, dimx_, dimf_;
  bool is_active_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_factorization.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_ 