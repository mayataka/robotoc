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

  ///
  /// @brief Set the dimension of the impulse. 
  /// @param[in] impulse_status Impulse status. 
  ///
  void setImpulseStatus(const ImpulseStatus& impulse_status);

  ///
  /// @brief A factorization matrix at a time stage. 
  /// @param[in] time_stage Time stage of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T(const int time_stage);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T(const int time_stage) const;

  ///
  /// @brief A factorization matrix at a impulse time stage. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_impulse(const int impulse_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T_impulse(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_impulse(
      const int impulse_index) const;

  ///
  /// @brief A factorization matrix at a aux time stage. 
  /// @param[in] impulse_index Impulse index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_aux(const int impulse_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T_aux(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_aux(const int impulse_index) const;

  ///
  /// @brief A factorization matrix at a lift time stage. 
  /// @param[in] lift_index Lift index of interested. 
  ///
  Eigen::Block<Eigen::MatrixXd> T_lift(const int lift_index);

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::T_lift(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> T_lift(const int lift_index) const;

  ///
  /// @brief Partial derivative of the equality constriant with respect to the 
  /// configuration. 
  ///
  Eigen::Block<Eigen::MatrixXd> Eq();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::Eq(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> Eq() const;

  ///
  /// @brief Product of the partial derivative of the equality constriant with 
  /// respect to state and RiccatiFactorization::N. 
  ///
  Eigen::Block<Eigen::MatrixXd> EN();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::EN(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> EN() const;

  ///
  /// @brief Top Robot::dimv() rows of RiccatiFactorization::EN(). 
  ///
  Eigen::Block<Eigen::MatrixXd> ENq();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::ENq(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> ENq() const;

  ///
  /// @brief Product of RiccatiFactorization::EN() and 
  /// StateConstraintRiccatiFactorization::Eq(). 
  ///
  Eigen::Block<Eigen::MatrixXd> ENEt();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::ENEt(). 
  ///
  const Eigen::Block<const Eigen::MatrixXd> ENEt() const;

  ///
  /// @brief Residual of the equality constriant.
  ///
  Eigen::VectorBlock<Eigen::VectorXd> e();

  ///
  /// @brief const version of StateConstraintRiccatiFactorization::e(). 
  ///
  const Eigen::VectorBlock<const Eigen::VectorXd> e() const;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

private:
  std::vector<Eigen::MatrixXd> T_full_, T_impulse_full_, T_aux_full_, T_lift_full_;
  Eigen::MatrixXd E_full_, EN_full_, ENEt_full_;
  Eigen::VectorXd e_full_;
  int N_, max_num_impulse_, dimv_, dimx_, dimf_;

};

} // namespace idocp

#include "idocp/ocp/state_constraint_riccati_factorization.hxx"

#endif // IDOCP_STATE_CONSTRAINT_RICCATI_FACTORIZATION_HPP_ 