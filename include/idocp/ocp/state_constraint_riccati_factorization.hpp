#ifndef IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_
#define IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_

#include <vector>

#include "Eigen/Core"

#include "idocp/robot/robot.hpp"
#include "idocp/robot/impulse_status.hpp"

namespace idocp {

///
/// @class StateConstraintRiccatiFactorization
/// @brief Riccati factorized matrix and vector for a state only equality 
/// constraint.
///
class StateConstraintRiccatiFactorization {
public:
  ///
  /// @brief Allocate Riccati factorization matrix and vector.
  /// @param[in] robot Robot model. Must be initialized by URDF or XML.
  /// @param[in] N Number of discretization of the horizon. Must be more than 1. 
  /// @param[in] max_num_impulse Maximum number of possible impulse. Must be 
  /// more than 1 and less than N. 
  ///
  StateConstraintRiccatiFactorization(const Robot& robot, const int N, 
                                      const int max_num_impulse);

  ///
  /// @brief Default constructor. 
  ///
  StateConstraintRiccatiFactorization();

  ///
  /// @brief Destructor. 
  ///
  ~StateConstraintRiccatiFactorization();
 
  ///
  /// @brief Default copy constructor. 
  ///
  StateConstraintRiccatiFactorization(
      const StateConstraintRiccatiFactorization&) = default;

  ///
  /// @brief Default copy operator. 
  ///
  StateConstraintRiccatiFactorization& operator=(
      const StateConstraintRiccatiFactorization&) = default;

  ///
  /// @brief Default move constructor. 
  ///
  StateConstraintRiccatiFactorization(
      StateConstraintRiccatiFactorization&&) noexcept = default;

  ///
  /// @brief Default move assign operator. 
  ///
  StateConstraintRiccatiFactorization& operator=(
      StateConstraintRiccatiFactorization&&) noexcept = default;

  void setImpulseStatus(const ImpulseStatus& impulse_status);

  Eigen::Block<Eigen::MatrixXd> T(const int time_stage);

  const Eigen::Block<const Eigen::MatrixXd> T(const int time_stage) const;

  Eigen::Block<Eigen::MatrixXd> T_impulse(const int impulse_index);

  const Eigen::Block<const Eigen::MatrixXd> T_impulse(const int impulse_index) const;

  Eigen::Block<Eigen::MatrixXd> T_lift(const int lift_index);

  const Eigen::Block<const Eigen::MatrixXd> T_lift(const int lift_index) const;

  Eigen::Block<Eigen::MatrixXd> Eq();

  const Eigen::Block<const Eigen::MatrixXd> Eq() const;

  Eigen::Block<Eigen::MatrixXd> EN();

  const Eigen::Block<const Eigen::MatrixXd> EN() const;

  Eigen::Block<Eigen::MatrixXd> ENq();

  const Eigen::Block<const Eigen::MatrixXd> ENq() const;

  Eigen::Block<Eigen::MatrixXd> ENEt();

  const Eigen::Block<const Eigen::MatrixXd> ENEt() const;

  Eigen::VectorBlock<Eigen::VectorXd> e();

  const Eigen::VectorBlock<const Eigen::VectorXd> e() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::MatrixXd> T_full_, T_impulse_full_, T_lift_full_;
  Eigen::MatrixXd E_full_, EN_full_, ENEt_full_;
  Eigen::VectorXd e_full_;
  int N_, max_num_impulse_, dimv_, dimx_, dimf_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_factorization.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_ 